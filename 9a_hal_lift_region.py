#!/usr/bin/env python3
"""
Lift the HelSus SLC4A8-AQP5 interval through the Cactus alignment into every other genome,
for the alignment-based synteny figure (9b_plot_hal_synteny.py).

1. Region. SLC4A8 and AQP5 are found in the HelSus annotation (--gff; EviAnn names such as
   'Similar to SLC4A8: ...' are read as SLC4A8; NAME=ID adds gene IDs, e.g.
   AQP5=LOC_00012326). The region is the pair +/- --pad-kb.
2. Tiles. The region is cut into --tile-kb tiles (BED6, + strand) and every HelSus gene in
   it is written as a BED6 interval.
3. Lift. halLiftover maps tiles and genes from HelSus into each target genome (every leaf
   of the HAL tree except HelSus, or --targets). A lifted '-' strand = inverted relative to
   HelSus.
4. Hits. Per tile and target: the contig with the most lifted bp, and on it the largest
   cluster of fragments (gap < --cluster-kb), is the primary hit. Other clusters with at
   least --secondary-frac of the tile's bp are secondary hits (duplicated or paralogous
   sequence).
5. Blocks. Consecutive hits on one contig and strand whose target positions continue in
   order (gap < --max-gap-kb in both genomes) are chained into a syntenic block.

Needs halLiftover and halStats on PATH (Biowulf: module load cactus). Python standard
library only. Read-only on the HAL file and GFF. Writes to --outdir:
  region.tsv         HelSus chromosome, region, boundary genes used
  ref_genes.tsv      HelSus genes in the region (symbol, description, family group)
  tiles.bed, genes.bed           what was lifted
  lifted/<target>.{tiles,genes}.bed   raw halLiftover output
  blocks.tsv         target, contig, target start/end, strand, HelSus start/end, aligned bp,
                     primary/secondary, number of tiles
  genes_lifted.tsv   target, HelSus gene, contig, start, end, strand, fraction of the gene lifted
  targets.tsv        target genome -> species code, aligned bp, main contig

Usage
  module load cactus
  python3 9a_hal_lift_region.py --hal ${CACTUSDIR}/<alignment>.hal --ref HelSus \
      --gff ${HELSUS_EVIANN_GFF} --boundaries SLC4A8,AQP5=LOC_00012326 \
      --pad-kb 25 --outdir hal_lift --threads 8 2> hal_lift.log
"""
import argparse, gzip, os, re, shutil, subprocess, sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import unquote

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

DEFAULT_RX = {"SLC4A8": r"solute carrier family 4 member 8", "AQP5": r"aquaporin[- ]?5"}
SIMILAR = re.compile(r"^\s*similar to\s+([A-Za-z0-9][A-Za-z0-9._-]*)\s*:\s*(.*)$", re.I)
ID_LIKE = re.compile(r"^(LOC_?\d+|gene-|rna-|[XN][MR]_\d+|ENS[A-Z]*\d+)", re.I)


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def code_of(name):
    """HAL genome name -> 3+3 code: 'HelSus' stays, 'Heloderma_suspectum' -> 'HelSus'."""
    if re.fullmatch(r"[A-Z][a-z]{2}[A-Z][a-z]{2}", name):
        return name
    parts = [p for p in re.split(r"[_\s.]+", name) if p]
    if len(parts) >= 2:
        return parts[0][:3].capitalize() + parts[1][:3].capitalize()
    return name


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"[error] {' '.join(cmd)}\n{r.stderr.strip()[:500]}")
    return r.stdout


def hal_leaves(hal):
    tree = run(["halStats", "--tree", hal]).strip()
    s = re.sub(r":[0-9.eE+-]+", "", tree)
    return [t.strip() for t in re.findall(r"[(,]\s*([^(),;]+)", s) if t.strip()]


def gff_genes(gff, contig=None):
    """gene records: dict(id, contig, start, end, strand, symbol, description)."""
    genes, txdesc = {}, {}
    with opener(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 9 or (contig and t[0] != contig):
                continue
            at = {k: unquote(v) for k, v in (kv.split("=", 1) for kv in t[8].strip().strip(";").split(";") if "=" in kv)}
            if t[2] in ("gene", "pseudogene"):
                sym = at.get("gene") or at.get("Name") or ""
                if sym == at.get("ID") or ID_LIKE.match(sym):
                    sym = ""
                desc = at.get("description") or at.get("Note") or at.get("product") or ""
                genes[at.get("ID", "")] = dict(id=at.get("ID", ""), contig=t[0], start=int(t[3]), end=int(t[4]),
                                               strand=t[6], symbol=sym, description=desc,
                                               biotype=at.get("gene_biotype", t[2]))
            elif t[2] in ("mRNA", "transcript") and at.get("Parent"):
                txdesc.setdefault(at["Parent"].split(",")[0], at.get("product") or at.get("Note") or "")
    for gid, g in genes.items():
        if not g["description"]:
            g["description"] = txdesc.get(gid, "")
        if not g["symbol"]:
            m = SIMILAR.match(g["description"])
            if m:
                g["symbol"], g["description"] = m.group(1), m.group(2)
        text = " | ".join([g["symbol"], g["description"]])
        grp = N.name_group(text)
        g["family"] = N.is_family(grp) or grp == "proglucagon"
        g["group"] = grp if g["family"] else "other gene"
    return list(genes.values())


def find_region(genes, bspec, pad):
    cands = {}
    for name, ids, rx in bspec:
        cands[name] = [g for g in genes if g["symbol"].upper() == name.upper() or g["id"] in ids
                       or rx.search(g["description"])]
        if not cands[name]:
            sys.exit(f"[error] no {name} gene in the GFF; give its ID as {name}=<gene ID>")
    (n0, *_), (n1, *_) = bspec
    pairs = [(abs((a["start"] + a["end"]) / 2 - (b["start"] + b["end"]) / 2), a, b)
             for a in cands[n0] for b in cands[n1] if a["contig"] == b["contig"]]
    if not pairs:
        sys.exit(f"[error] {n0} and {n1} are not on one contig in the GFF")
    _, a, b = min(pairs, key=lambda p: p[0])
    for name, c in ((n0, cands[n0]), (n1, cands[n1])):
        if len(c) > 1:
            N.info(f"{name}: {len(c)} candidates ({', '.join(g['id'] for g in c[:6])}); using the pair closest together")
    lo, hi = min(a["start"], b["start"]), max(a["end"], b["end"])
    return a["contig"], max(1, lo - pad), hi + pad, a, b


def read_bed(path):
    out = defaultdict(list)
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) >= 4:
                out[t[3]].append((t[0], int(t[1]), int(t[2]), t[5] if len(t) >= 6 else "+"))
    return out


def clusters(frags, gap):
    frags = sorted(frags, key=lambda f: f[1])
    out, cur = [], []
    for f in frags:
        if cur and f[1] - max(x[2] for x in cur) > gap:
            out.append(cur); cur = []
        cur.append(f)
    if cur:
        out.append(cur)
    return out


def summarise(frags, gap):
    """fragments of one item -> [(contig, start, end, strand, bp)] best first."""
    by_ctg = defaultdict(list)
    for f in frags:
        by_ctg[f[0]].append(f)
    hits = []
    for c, fs in by_ctg.items():
        for cl in clusters(fs, gap):
            bp = sum(f[2] - f[1] for f in cl)
            plus = sum(f[2] - f[1] for f in cl if f[3] == "+")
            hits.append((c, min(f[1] for f in cl), max(f[2] for f in cl), "+" if plus * 2 >= bp else "-", bp))
    return sorted(hits, key=lambda h: -h[4])


def chain(hits, max_gap):
    """hits: [(rstart, rend, contig, tstart, tend, strand, bp)] in HelSus order -> blocks."""
    blocks = []
    for h in sorted(hits):
        rs, re_, c, ts, te, st, bp = h
        b = blocks[-1] if blocks else None
        ok = False
        if b and b["contig"] == c and b["strand"] == st and rs - b["rend"] <= max_gap:
            if st == "+":
                ok = -max_gap / 4 <= ts - b["tend"] <= max_gap
            else:
                ok = -max_gap / 4 <= b["tstart"] - te <= max_gap
        if ok:
            b["rend"] = max(b["rend"], re_); b["tstart"] = min(b["tstart"], ts); b["tend"] = max(b["tend"], te)
            b["bp"] += bp; b["n"] += 1
        else:
            blocks.append(dict(contig=c, tstart=ts, tend=te, strand=st, rstart=rs, rend=re_, bp=bp, n=1))
    return blocks


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--hal", required=True)
    ap.add_argument("--ref", default="HelSus", help="reference genome: HAL name, or its 3+3 code")
    ap.add_argument("--gff", required=True, help="annotation of the reference (HelSus EviAnn)")
    ap.add_argument("--hal-chrom", help="reference chromosome name in the HAL, if not the GFF's")
    ap.add_argument("--boundaries", default="SLC4A8,AQP5", help="left,right; NAME=ID|ID adds gene IDs")
    ap.add_argument("--pad-kb", type=float, default=25.0, help="added to each side of the boundary pair")
    ap.add_argument("--tile-kb", type=float, default=1.0)
    ap.add_argument("--targets", help="comma-separated HAL genome names (default: every other leaf)")
    ap.add_argument("--cluster-kb", type=float, default=5.0, help="fragments of one tile/gene closer than this are one hit")
    ap.add_argument("--secondary-frac", type=float, default=0.3, help="secondary hits need this fraction of the primary's bp")
    ap.add_argument("--max-gap-kb", type=float, default=50.0, help="largest gap inside a syntenic block")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    for tool in ("halLiftover", "halStats"):
        if not shutil.which(tool):
            sys.exit(f"[error] {tool} not on PATH (module load cactus)")
    os.makedirs(os.path.join(a.outdir, "lifted"), exist_ok=True)

    leaves = hal_leaves(a.hal)
    ref = a.ref if a.ref in leaves else next((l for l in leaves if code_of(l) == a.ref), None)
    if not ref:
        sys.exit(f"[error] reference {a.ref} is not a HAL leaf: {', '.join(leaves)}")
    targets = a.targets.split(",") if a.targets else [l for l in leaves if l != ref]
    N.info(f"HAL leaves: {len(leaves)}; reference {ref}; {len(targets)} targets")

    bspec = []
    for item in a.boundaries.split(","):
        name, _, ids = item.partition("=")
        bspec.append((name, {i for i in ids.split("|") if i},
                      re.compile(rf"(?<![A-Za-z0-9])(?:{DEFAULT_RX.get(name, re.escape(name))})(?![A-Za-z0-9])", re.I)))
    genes = gff_genes(a.gff)
    contig, lo, hi, b0, b1 = find_region(genes, bspec, int(a.pad_kb * 1000))
    chrom = a.hal_chrom or contig
    seqs = set(run(["halStats", "--sequences", ref, a.hal]).replace(",", "\n").split())
    if chrom not in seqs:
        sys.exit(f"[error] {chrom} is not a sequence of {ref} in the HAL (e.g. {', '.join(sorted(seqs)[:5])}); "
                 f"pass the HAL name with --hal-chrom")
    N.info(f"region {chrom}:{lo}-{hi} ({(hi - lo) / 1000:.0f} kb): {bspec[0][0]} {b0['id']} {b0['start']}-{b0['end']}, "
           f"{bspec[1][0]} {b1['id']} {b1['start']}-{b1['end']}, +/- {a.pad_kb:g} kb")
    with open(os.path.join(a.outdir, "region.tsv"), "w") as fh:
        fh.write("ref\tchrom\tstart\tend\tleft\tleft_id\tright\tright_id\n")
        fh.write(f"{ref}\t{chrom}\t{lo}\t{hi}\t{bspec[0][0]}\t{b0['id']}\t{bspec[1][0]}\t{b1['id']}\n")

    rgenes = [g for g in genes if g["contig"] == contig and g["end"] >= lo and g["start"] <= hi]
    for g in rgenes:
        g["kind"] = "boundary" if g["id"] in (b0["id"], b1["id"]) else "gene"
    with open(os.path.join(a.outdir, "ref_genes.tsv"), "w") as fh:
        fh.write("gene_id\tsymbol\tgroup\tfamily\tkind\tstart\tend\tstrand\tdescription\n")
        for g in sorted(rgenes, key=lambda g: g["start"]):
            fh.write(f"{g['id']}\t{g['symbol']}\t{g['group']}\t{int(g['family'])}\t{g['kind']}\t{g['start']}\t"
                     f"{g['end']}\t{g['strand']}\t{g['description']}\n")
    N.info(f"{len(rgenes)} HelSus genes in the region, {sum(g['family'] for g in rgenes)} exendin-family")

    tile = int(a.tile_kb * 1000)
    tiles_bed, genes_bed = os.path.join(a.outdir, "tiles.bed"), os.path.join(a.outdir, "genes.bed")
    with open(tiles_bed, "w") as fh:
        for s in range(lo - 1, hi, tile):
            fh.write(f"{chrom}\t{s}\t{min(s + tile, hi)}\tt{s}\t0\t+\n")
    with open(genes_bed, "w") as fh:
        for g in rgenes:
            fh.write(f"{chrom}\t{g['start'] - 1}\t{g['end']}\t{g['id']}\t0\t+\n")

    def lift(tg):
        outs = []
        for kind, src in (("tiles", tiles_bed), ("genes", genes_bed)):
            dst = os.path.join(a.outdir, "lifted", f"{tg}.{kind}.bed")
            r = subprocess.run(["halLiftover", a.hal, ref, src, tg, dst], capture_output=True, text=True)
            outs.append(r.returncode == 0 or print(f"[warn] halLiftover {tg} {kind}: {r.stderr.strip()[:300]}",
                                                   file=sys.stderr))
        return tg, all(outs)

    with ThreadPoolExecutor(max_workers=a.threads) as ex:
        for tg, ok in ex.map(lift, targets):
            N.info(f"lifted to {tg}" + ("" if ok else " (with errors)"))

    gap, maxgap = a.cluster_kb * 1000, a.max_gap_kb * 1000
    glen = {g["id"]: g["end"] - g["start"] + 1 for g in rgenes}
    with open(os.path.join(a.outdir, "blocks.tsv"), "w") as fb, \
         open(os.path.join(a.outdir, "genes_lifted.tsv"), "w") as fg, \
         open(os.path.join(a.outdir, "targets.tsv"), "w") as ft:
        fb.write("target\tcode\tcontig\ttstart\ttend\tstrand\trstart\trend\taligned_bp\ttype\tn_tiles\n")
        fg.write("target\tcode\tgene_id\tcontig\tstart\tend\tstrand\tfrac_lifted\thit\n")
        ft.write("target\tcode\taligned_bp\tfrac_region\tmain_contig\tn_blocks\n")
        for tg in targets:
            code = code_of(tg)
            lifted = read_bed(os.path.join(a.outdir, "lifted", f"{tg}.tiles.bed"))
            prim, sec = [], []
            for name, frags in lifted.items():
                rs = int(name[1:]) + 1
                re_ = min(rs + tile - 1, hi)
                hs = [h for h in summarise(frags, gap) if h[4] > 0]
                if not hs:
                    continue
                c, s, e, st, bp = hs[0]
                prim.append((rs, re_, c, s, e, st, bp))
                for c2, s2, e2, st2, bp2 in hs[1:]:
                    if bp2 >= a.secondary_frac * bp:
                        sec.append((rs, re_, c2, s2, e2, st2, bp2))
            pb, sb = chain(prim, maxgap), chain(sec, maxgap)
            for kind, bl in (("primary", pb), ("secondary", sb)):
                for b in bl:
                    fb.write(f"{tg}\t{code}\t{b['contig']}\t{b['tstart']}\t{b['tend']}\t{b['strand']}\t{b['rstart']}\t"
                             f"{b['rend']}\t{b['bp']}\t{kind}\t{b['n']}\n")
            bp_ctg = defaultdict(int)
            for b in pb:
                bp_ctg[b["contig"]] += b["bp"]
            tot = sum(bp_ctg.values())
            main = max(bp_ctg, key=bp_ctg.get) if bp_ctg else "NA"
            ft.write(f"{tg}\t{code}\t{tot}\t{tot / (hi - lo + 1):.3f}\t{main}\t{len(pb)}\n")
            N.info(f"{tg} ({code}): {tot / 1000:.0f} kb of {(hi - lo) / 1000:.0f} kb aligned "
                   f"({100 * tot / (hi - lo + 1):.0f}%), {len(pb)} primary block(s) on {len(bp_ctg)} contig(s), "
                   f"main {main}; {len(sb)} secondary block(s)")
            glift = read_bed(os.path.join(a.outdir, "lifted", f"{tg}.genes.bed"))
            for gid, frags in glift.items():
                hs = summarise(frags, max(gap, 20000))
                for k, (c, s, e, st, bp) in enumerate(hs[:3]):
                    if k and bp < a.secondary_frac * hs[0][4]:
                        break
                    fg.write(f"{tg}\t{code}\t{gid}\t{c}\t{s + 1}\t{e}\t{st}\t{bp / glen.get(gid, 1):.3f}\t"
                             f"{'primary' if k == 0 else 'secondary'}\n")
    N.info(f"wrote {a.outdir}/blocks.tsv, genes_lifted.tsv, targets.tsv")


if __name__ == "__main__":
    main()
