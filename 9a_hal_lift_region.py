#!/usr/bin/env python3
"""
Lift the HelSus SLC4A8-AQP5 interval through the Cactus alignment into every other genome,
for the alignment-based river plot (9b_plot_hal_river.py).

Genes come from the TOGA2 interval tables written by 3c (--genes-out), not from a GFF, so
every genome carries ITS OWN annotation in its own coordinates rather than a projection of
the reference's. --names applies the step-5 nomenclature (Exendin_anc, A, B, Z, GCGL ...).

1. Region. The reference's two boundary rows (the `boundary` column of the TOGA2 table)
   give SLC4A8 and AQP5; the region is that pair +/- --pad-kb.
2. Tiles. The region is cut into --tile-kb tiles (BED6, + strand). The reference's own genes
   are also written as BED6, so their projection can be compared with each target's
   annotation.
3. Lift. halLiftover maps tiles and reference genes from the reference into each target
   genome (every leaf of the HAL tree except the reference, or --targets). A lifted '-'
   strand = inverted relative to the reference.
4. Hits. Per tile and target: the contig with the most lifted bp, and on it the largest
   cluster of fragments (gap < --cluster-kb), is the primary hit. Other clusters with at
   least --secondary-frac of the tile's bp are secondary hits (duplicated or paralogous
   sequence).
5. Blocks. Consecutive hits on one contig and strand whose target positions continue in
   order (gap < --max-gap-kb in both genomes) are chained into a syntenic block.

Needs halLiftover and halStats on PATH (Biowulf: module load cactus). Python standard
library only. Read-only on the HAL file and the TOGA2 tables. Writes to --outdir:
  region.tsv         reference chromosome, region, boundary genes used
  ref_genes.tsv      reference genes in the region (name, symbol, group, family, kind)
  target_genes.tsv   every target's OWN TOGA2 genes, in its own coordinates, with the
                     step-5 name where --names gives one
  tiles.bed, genes.bed           what was lifted
  lifted/<target>.{tiles,genes}.bed   raw halLiftover output
  tile_hits.tsv      every tile's primary and secondary hit in every target (for the river plot)
  blocks.tsv         target, contig, target start/end, strand, reference start/end, aligned bp,
                     primary/secondary, number of tiles
  genes_lifted.tsv   target, reference gene, contig, start, end, strand, fraction lifted
  targets.tsv        target genome -> species code, aligned bp, main contig

Usage
  module load cactus
  python3 9a_hal_lift_region.py --hal ${CACTUSDIR}/<alignment>.hal --ref HelSus \
      --toga2-genes ../synteny_toga2/interval/interval_HelSus_genes.tsv \
      --names ../toga2_exendin_tree/names/names.tsv \
      --pad-kb 25 --outdir hal_lift --threads 8 2> hal_lift.log
"""
import argparse, csv, os, re, shutil, subprocess, sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import unquote

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

def read_tsv(path):
    if not os.path.exists(path):
        sys.exit(f"[error] no such file: {path}")
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


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


def toga2_genes(paths, names_path=None):
    """{species code: [gene records]} from the 3c --genes-out tables. First file wins per species."""
    byspp, src = {}, {}
    for p in paths:
        for r in read_tsv(p):
            q = r.get("species", "")
            if q in src and src[q] != p:
                continue
            src.setdefault(q, p)
            grp = r.get("group") or ""
            fam = bool(N.is_family(grp) or grp == "proglucagon")
            byspp.setdefault(q, []).append(dict(
                id=r.get("query_gene", ""), contig=r.get("contig", ""),
                start=int(r["start"]), end=int(r["end"]), strand=r.get("strand") or "+",
                symbol=r.get("symbol") or "", description=r.get("description") or "",
                status=r.get("status") or "", group=grp if fam else "other gene", family=fam,
                kind="boundary" if r.get("boundary") else "gene", name=""))
    for q in sorted(src):
        N.info(f"{q}: {len(byspp[q])} TOGA2 genes ({sum(g['family'] for g in byspp[q])} family) "
               f"from {os.path.basename(src[q])}")
    if names_path:
        apply_names(byspp, names_path)
    return byspp


def apply_names(byspp, path):
    """names.tsv from 10f/10h -> the 'name' field, matched on species + contig + overlap."""
    n, miss = 0, []
    for r in read_tsv(path):
        q, new = r.get("species", ""), r.get("new_name", "")
        try:
            s0, e0 = int(r["start"]), int(r["end"])
        except (KeyError, ValueError):
            miss.append(new); continue
        best = None
        for g in byspp.get(q, []):
            if g["contig"] != r.get("contig"):
                continue
            ov = min(g["end"], e0) - max(g["start"], s0)
            if ov > 0 and (best is None or ov > best[0]):
                best = (ov, g)
        if best:
            best[1]["name"] = new; n += 1
        else:
            miss.append(new)
    N.info(f"named {n} gene(s) from {os.path.basename(path)}")
    if miss:
        N.warn(f"{len(miss)} name(s) with no gene in the TOGA2 interval table: {', '.join(miss[:8])}"
               + (" ..." if len(miss) > 8 else ""))


def find_region(genes, bnames, pad):
    """the two boundary rows of the reference -> contig, lo, hi, left gene, right gene."""
    cands = {}
    for nm in bnames:
        cands[nm] = [g for g in genes if g["kind"] == "boundary" and g["symbol"].upper() == nm.upper()]
        if not cands[nm]:
            have = sorted({g["symbol"] for g in genes if g["kind"] == "boundary"})
            sys.exit(f"[error] no {nm} boundary row for the reference in the TOGA2 table "
                     f"(boundary symbols present: {', '.join(have) or 'none'}); check --boundaries")
    n0, n1 = bnames
    pairs = [(abs((x["start"] + x["end"]) / 2 - (y["start"] + y["end"]) / 2), x, y)
             for x in cands[n0] for y in cands[n1] if x["contig"] == y["contig"]]
    if not pairs:
        sys.exit(f"[error] {n0} and {n1} are not on one contig in the reference's TOGA2 rows")
    _, x, y = min(pairs, key=lambda t: t[0])
    for nm in bnames:
        if len(cands[nm]) > 1:
            N.info(f"{nm}: {len(cands[nm])} boundary rows; using the pair closest together")
    lo, hi = min(x["start"], y["start"]), max(x["end"], y["end"])
    return x["contig"], max(1, lo - pad), hi + pad, x, y


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
    ap.add_argument("--toga2-genes", required=True,
                    help="comma-separated 3c --genes-out tables; the first one holding a species supplies it")
    ap.add_argument("--names", help="names.tsv from 10f/10h: applies the step-5 gene names")
    ap.add_argument("--hal-chrom", help="reference chromosome name in the HAL, if not the table's")
    ap.add_argument("--boundaries", default="SLC4A8,AQP5", help="left,right boundary gene symbols")
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

    bnames = [x.strip() for x in a.boundaries.split(",")][:2]
    byspp = toga2_genes([x.strip() for x in a.toga2_genes.split(",") if x.strip()], a.names)
    refcode = code_of(ref)
    if refcode not in byspp:
        sys.exit(f"[error] the reference {refcode} has no rows in the TOGA2 table(s) "
                 f"(species present: {', '.join(sorted(byspp))})")
    genes = byspp[refcode]
    contig, lo, hi, b0, b1 = find_region(genes, bnames, int(a.pad_kb * 1000))
    chrom = a.hal_chrom or contig
    seqs = set(run(["halStats", "--sequences", ref, a.hal]).replace(",", "\n").split())
    if chrom not in seqs:
        sys.exit(f"[error] {chrom} is not a sequence of {ref} in the HAL (e.g. {', '.join(sorted(seqs)[:5])}); "
                 f"pass the HAL name with --hal-chrom")
    N.info(f"region {chrom}:{lo}-{hi} ({(hi - lo) / 1000:.0f} kb): {bnames[0]} {b0['id']} {b0['start']}-{b0['end']}, "
           f"{bnames[1]} {b1['id']} {b1['start']}-{b1['end']}, +/- {a.pad_kb:g} kb")
    with open(os.path.join(a.outdir, "region.tsv"), "w") as fh:
        fh.write("ref\tchrom\tstart\tend\tleft\tleft_id\tright\tright_id\n")
        fh.write(f"{ref}\t{chrom}\t{lo}\t{hi}\t{bnames[0]}\t{b0['id']}\t{bnames[1]}\t{b1['id']}\n")

    rgenes = [g for g in genes if g["contig"] == contig and g["end"] >= lo and g["start"] <= hi]
    for g in rgenes:
        g["kind"] = "boundary" if g["id"] in (b0["id"], b1["id"]) else "gene"
    with open(os.path.join(a.outdir, "ref_genes.tsv"), "w") as fh:
        fh.write("gene_id\tname\tsymbol\tgroup\tfamily\tkind\tstart\tend\tstrand\tstatus\tdescription\n")
        for g in sorted(rgenes, key=lambda g: g["start"]):
            fh.write(f"{g['id']}\t{g['name']}\t{g['symbol']}\t{g['group']}\t{int(g['family'])}\t{g['kind']}\t"
                     f"{g['start']}\t{g['end']}\t{g['strand']}\t{g['status']}\t{g['description']}\n")
    N.info(f"{len(rgenes)} {refcode} genes in the region, {sum(g['family'] for g in rgenes)} exendin-family")

    # every target's own annotation, in its own coordinates
    tcodes = {code_of(t): t for t in targets}
    with open(os.path.join(a.outdir, "target_genes.tsv"), "w") as fh:
        fh.write("code\tgene_id\tname\tsymbol\tgroup\tfamily\tkind\tcontig\tstart\tend\tstrand\tstatus\tdescription\n")
        nt = 0
        for q in sorted(tcodes):
            for g in sorted(byspp.get(q, []), key=lambda g: (g["contig"], g["start"])):
                fh.write(f"{q}\t{g['id']}\t{g['name']}\t{g['symbol']}\t{g['group']}\t{int(g['family'])}\t"
                         f"{g['kind']}\t{g['contig']}\t{g['start']}\t{g['end']}\t{g['strand']}\t{g['status']}\t"
                         f"{g['description']}\n")
                nt += 1
    no_ann = sorted(q for q in tcodes if not byspp.get(q))
    N.info(f"{nt} TOGA2 genes across {len(tcodes) - len(no_ann)} target(s)")
    if no_ann:
        N.warn(f"no TOGA2 rows (the river will show these rows without genes): {', '.join(no_ann)}")

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
        fh_t = open(os.path.join(a.outdir, "tile_hits.tsv"), "w")
        fh_t.write("target\tcode\trstart\trend\ttype\tcontig\ttstart\ttend\tstrand\taligned_bp\n")
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
            for kind, hl in (("primary", prim), ("secondary", sec)):
                for rs, re_, c, ts, te, st, bp in sorted(hl):
                    fh_t.write(f"{tg}\t{code}\t{rs}\t{re_}\t{kind}\t{c}\t{ts + 1}\t{te}\t{st}\t{bp}\n")
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
        fh_t.close()
    N.info(f"wrote {a.outdir}/blocks.tsv, tile_hits.tsv, genes_lifted.tsv, target_genes.tsv, targets.tsv")


if __name__ == "__main__":
    main()
