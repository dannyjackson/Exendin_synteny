#!/usr/bin/env python3
"""
Step 13: halSynteny PSL -> GENESPACE block table, for riparian plots of the Cactus alignment.

GENESPACE has no importer for whole-genome-alignment synteny, but its plotting core,
riparian_engine(), draws any block table with the columns genome1, genome2, chr1, chr2,
startBp1, endBp1, startBp2, endBp2, blkID (and a color). This script only reformats
halSynteny's PSL into that table; the blocks themselves are halSynteny's.

Two subcommands.

  genomes   Row order and HAL names for the figure, from the step-9 lift (9a --outdir) and
            the species tree, in the same order as the step-9 river plot. Writes a TSV:
              rank  code  hal_genome  interval_contigs
            interval_contigs = the HAL contigs the HelSus SLC4A8-AQP5 interval aligns to in
            that genome (9a blocks.tsv), for the fast --queryChromosome halSynteny runs.

  blk       halSynteny PSLs for each neighbouring pair of rows -> GENESPACE blk CSV.
            PSL files are <psl-dir>/<halA>__<halB>*.psl, run with --queryGenome halA
            --targetGenome halB, where halA is the upper row. Several files per pair
            (one per --queryChromosome) are concatenated.
            * refChr: each block is labelled with the HelSus chromosome it traces back to,
              chained outward from HelSus through the neighbouring pairs, so braids can be
              coloured by HelSus chromosome as in GENESPACE.
            * highlight: the HelSus interval (9a region.tsv) is projected through the same
              chain, block by block (linear within a block), and written as extra rows with
              highlight=1, drawn on top in the highlight colour.
            * --window-kb: zoom mode. Only the projected interval +/- this much is kept in
              every genome; blocks are clipped to those windows and coordinates are shifted
              so each window starts at 1. Contigs are renamed <contig>:<start Mb>.
            Inverted blocks (strand +- / -+) are written with startBp2 > endBp2, which is how
            riparian_engine draws the twist.

Standard library only. Read-only on its inputs.

Usage
  python3 13a_hal_to_blk.py genomes --indir ../synteny_cactus/hal_lift \
      --species-tree ${SPTREE} --species-table ../species_toga2.tsv --only-listed \
      --drop-species AnnSte,BarImb,ElgCoe --out genomes.tsv
  python3 13a_hal_to_blk.py blk --genomes genomes.tsv --psl-dir psl_genome \
      --region ../synteny_cactus/hal_lift/region.tsv --min-chr-mb 1 --out blk_genome.csv
  python3 13a_hal_to_blk.py blk --genomes genomes.tsv --psl-dir psl_region \
      --region ../synteny_cactus/hal_lift/region.tsv --window-kb 2000 --out blk_region.csv
"""
import argparse, bisect, csv, glob, os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def info(m):
    print(f"[info] {m}", file=sys.stderr)


def warn(m):
    print(f"[warn] {m}", file=sys.stderr)


def read_tsv(path):
    if not os.path.exists(path):
        sys.exit(f"[error] no such file: {path}")
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def code_of(name):
    """HAL genome name -> 3+3 code (as 9a): 'HelSus' stays, 'Heloderma_suspectum' -> 'HelSus'."""
    if re.fullmatch(r"[A-Z][a-z]{2}[A-Z][a-z]{2}", name):
        return name
    parts = [p for p in re.split(r"[_\s.]+", name) if p]
    return parts[0][:3].capitalize() + parts[1][:3].capitalize() if len(parts) >= 2 else name


# ------------------------------------------------------------------ genomes
def cmd_genomes(a):
    reg = read_tsv(os.path.join(a.indir, "region.tsv"))[0]
    refhal, refcode = reg["ref"], code_of(reg["ref"])
    hal = {refcode: refhal}
    for t in read_tsv(os.path.join(a.indir, "targets.tsv")):
        hal[t["code"]] = t["target"]
    ctgs = defaultdict(list)
    ctgs[refcode] = [reg["chrom"]]
    for b in read_tsv(os.path.join(a.indir, "blocks.tsv")):
        if b["type"] == "primary" and b["contig"] not in ctgs[b["code"]]:
            ctgs[b["code"]].append(b["contig"])
    codes = list(hal)
    drop = {c.strip() for c in a.drop_species.split(",") if c.strip()}
    if a.only_listed:
        import lib_exendin_names as N
        rows = N.read_species_table(a.species_table)
        tq = [r for r in rows if (r.get("toga2_query") or "").strip()]
        listed = {r["code"] for r in (tq or rows)} | {refcode}
        off = sorted(c for c in codes if c not in listed)
        if off:
            info(f"not in the species table, left out: {', '.join(off)}")
        codes = [c for c in codes if c in listed]
    codes = [c for c in codes if c not in drop or c == refcode]
    if a.species_tree:
        import lib_exendin_names as N
        S = N.load_engine()
        order, _ = S.species_order(codes, refcode, a.species_tree, a.tree_map, None)
    else:
        order = [refcode] + [c for c in codes if c != refcode]
    with open(a.out, "w") as fh:
        fh.write("rank\tcode\thal_genome\tinterval_contigs\n")
        for i, c in enumerate(order, 1):
            fh.write(f"{i}\t{c}\t{hal[c]}\t{','.join(ctgs.get(c, []))}\n")
    nint = sum(1 for c in order if ctgs.get(c))
    info(f"{len(order)} genomes in tree order, {len(order) - 1} neighbouring pairs; "
         f"{nint} with the interval aligned (step 9); reference {refcode} = {refhal}")
    no = [c for c in order if not ctgs.get(c)]
    if no:
        warn(f"interval not aligned in step 9 (empty in the zoom plot): {', '.join(no)}")


# ------------------------------------------------------------------ blk
class Index:
    """intervals per chromosome, for overlap queries."""

    def __init__(self, items):            # items: (chr, s, e, payload)
        self.by = defaultdict(list)
        for c, s, e, p in items:
            self.by[c].append((s, e, p))
        self.starts, self.maxlen = {}, {}
        for c, v in self.by.items():
            v.sort(key=lambda x: x[0])
            self.starts[c] = [x[0] for x in v]
            self.maxlen[c] = max(x[1] - x[0] for x in v) + 1

    def overlaps(self, c, s, e):
        v = self.by.get(c)
        if not v:
            return []
        i = bisect.bisect_right(self.starts[c], e)
        out, lo = [], s - self.maxlen[c]
        while i > 0:
            i -= 1
            x = v[i]
            if x[0] < lo:
                break
            if x[1] >= s:
                out.append(x)
        return out


def read_psl(path):
    out = []
    with open(path) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) < 17 or not t[0].isdigit():
                continue                       # psLayout header or blank
            st = t[8]
            orient = "+" if st.count("-") % 2 == 0 else "-"     # ++ / -- same, +- / -+ inverted
            out.append(dict(c1=t[9], s1=int(t[11]) + 1, e1=int(t[12]), size1=int(t[10]),
                            c2=t[13], s2=int(t[15]) + 1, e2=int(t[16]), size2=int(t[14]),
                            o=orient, match=int(t[0])))
    return out


def project(b, s, e, side):
    """clip block b to [s, e] on `side` (1 or 2) -> (s1, e1, s2, e2) or None."""
    k, o_ = ("1", "2") if side == 1 else ("2", "1")
    bs, be = b["s" + k], b["e" + k]
    cs, ce = max(bs, s), min(be, e)
    if ce < cs:
        return None
    os_, oe = b["s" + o_], b["e" + o_]
    L, M = max(be - bs, 1), oe - os_

    def f(x):
        r = (x - bs) / L
        return os_ + r * M if b["o"] == "+" else oe - r * M
    p, q = sorted((round(f(cs)), round(f(ce))))
    return (cs, ce, p, q) if side == 1 else (p, q, cs, ce)


def merge(pieces, gap):
    """[(chr, s, e)] -> merged per chromosome where closer than gap."""
    out = []
    for c, s, e in sorted(pieces):
        if out and out[-1][0] == c and s - out[-1][2] <= gap:
            out[-1] = (c, out[-1][1], max(out[-1][2], e))
        else:
            out.append((c, s, e))
    return out


def cmd_blk(a):
    g = read_tsv(a.genomes)
    order = [r["code"] for r in g]
    halof = {r["code"]: r["hal_genome"] for r in g}
    reg = read_tsv(a.region)[0] if a.region else None
    ref = a.ref or (code_of(reg["ref"]) if reg else order[0])
    if ref not in order:
        sys.exit(f"[error] reference {ref} is not in {a.genomes}")
    ri = order.index(ref)

    # one block list per neighbouring pair, genome1 = upper row = halSynteny query
    pairs = []
    for up, dn in zip(order, order[1:]):
        files = sorted(glob.glob(os.path.join(a.psl_dir, f"{halof[up]}__{halof[dn]}*.psl")))
        bl = [b for f in files for b in read_psl(f)]
        if a.min_block_kb:
            bl = [b for b in bl if b["e1"] - b["s1"] + 1 >= a.min_block_kb * 1000]
        if not files:
            warn(f"{up} -> {dn}: no PSL {halof[up]}__{halof[dn]}*.psl in {a.psl_dir}; the river breaks here")
        pairs.append(dict(g1=up, g2=dn, blocks=bl))

    # chain HelSus chromosome labels and the interval outward from the reference
    lab = {ref: None}                           # genome -> Index of (chr, s, e, refChr); None = itself
    proj = {ref: [(reg["chrom"], int(reg["start"]), int(reg["end"]))] if reg else []}
    hl_rows = defaultdict(list)                 # pair index -> clipped interval rows

    def label_of(gid, c, s, e):
        if lab.get(gid) is None:
            return c if gid == ref else ""
        bp = defaultdict(int)
        for xs, xe, rc in lab[gid].overlaps(c, s, e):
            bp[rc] += min(xe, e) - max(xs, s) + 1
        return max(bp, key=bp.get) if bp else ""

    walk = [(i, 1) for i in range(ri, len(pairs))] + [(i, 2) for i in range(ri - 1, -1, -1)]
    for i, known in walk:                       # known side: 1 going down from ref, 2 going up
        p = pairs[i]
        kg, og = (p["g1"], p["g2"]) if known == 1 else (p["g2"], p["g1"])
        k, o = ("1", "2") if known == 1 else ("2", "1")
        items = []
        for b in p["blocks"]:
            b["refChr"] = label_of(kg, b["c" + k], b["s" + k], b["e" + k])
            if b["refChr"]:
                items.append((b["c" + o], b["s" + o], b["e" + o], b["refChr"]))
        lab[og] = Index(items)
        newp = []
        pidx = Index([(c, s, e, None) for c, s, e in proj.get(kg, [])])
        for b in p["blocks"]:
            for s, e, _ in pidx.overlaps(b["c" + k], b["s" + k], b["e" + k]):
                r = project(b, s, e, known)
                if r:
                    hl_rows[i].append((b, r))
                    newp.append((b["c" + o], r[0] if known == 2 else r[2], r[1] if known == 2 else r[3]))
        proj[og] = merge(newp, a.merge_kb * 1000)

    for c in order:
        pc = proj.get(c, [])
        if reg:
            info(f"{c}: interval projected to " + (", ".join(f"{x}:{s}-{e}" for x, s, e in pc[:4])
                                                   + (" ..." if len(pc) > 4 else "") if pc else "nothing"))

    # zoom: windows around the projected interval, blocks clipped to them
    rename = {}
    if a.window_kb is not None:
        if not reg:
            sys.exit("[error] --window-kb needs --region")
        w = int(a.window_kb * 1000)
        win = {}
        for c in order:
            ws = merge([(x, max(1, s - w), e + w) for x, s, e in proj.get(c, [])], 0)
            win[c] = Index([(x, s, e, (x, s, e)) for x, s, e in ws])
            for x, s, e in ws:
                rename[(c, x, s)] = f"{x}:{s / 1e6:.1f}"

        def clip(b, g1, g2):
            out = []
            for s, e, w1 in win[g1].overlaps(b["c1"], b["s1"], b["e1"]):
                r = project(b, s, e, 1)
                if not r:
                    continue
                sub = dict(b, s1=r[0], e1=r[1], s2=r[2], e2=r[3])
                for s_, e_, w2 in win[g2].overlaps(b["c2"], r[2], r[3]):
                    r2 = project(sub, s_, e_, 2)
                    if r2:
                        out.append((dict(b, s1=r2[0], e1=r2[1], s2=r2[2], e2=r2[3]), w1, w2))
            return out
    rows, uid = [], [0]

    def emit(p, b, w1=None, w2=None, hl=0):
        uid[0] += 1
        c1, c2 = b["c1"], b["c2"]
        s1, e1, s2, e2 = b["s1"], b["e1"], b["s2"], b["e2"]
        if w1:
            c1, s1, e1 = rename[(p["g1"], w1[0], w1[1])], s1 - w1[1] + 1, e1 - w1[1] + 1
            c2, s2, e2 = rename[(p["g2"], w2[0], w2[1])], s2 - w2[1] + 1, e2 - w2[1] + 1
        if b["o"] == "-":
            s2, e2 = e2, s2                     # riparian_engine draws start2 > end2 as a twist
        rows.append([p["g1"], p["g2"], c1, c2, s1, e1, s2, e2, b["o"],
                     f"{'hl' if hl else 'b'}{uid[0]}_{p['g1']}_{p['g2']}", ref, b.get("refChr", ""), hl])

    for i, p in enumerate(pairs):
        todo = [(b, 0) for b in p["blocks"]] + \
               [(dict(b, s1=r[0], e1=r[1], s2=r[2], e2=r[3]), 1) for b, r in hl_rows[i]]
        for b, hl in todo:
            if a.window_kb is None:
                emit(p, b, hl=hl)
            else:
                for cb, w1, w2 in clip(b, p["g1"], p["g2"]):
                    emit(p, cb, w1, w2, hl=hl)

    # drop chromosomes shorter than --min-chr-mb (as drawn: up to the last block)
    if a.min_chr_mb:
        ext = defaultdict(int)
        for r in rows:
            ext[(r[0], r[2])] = max(ext[(r[0], r[2])], r[4], r[5])
            ext[(r[1], r[3])] = max(ext[(r[1], r[3])], r[6], r[7])
        keep = {k for k, v in ext.items() if v >= a.min_chr_mb * 1e6}
        before = len(rows)
        rows = [r for r in rows if (r[0], r[2]) in keep and (r[1], r[3]) in keep]
        info(f"--min-chr-mb {a.min_chr_mb:g}: {len(ext) - len(keep)} contig(s) and "
             f"{before - len(rows)} block(s) left out")

    with open(a.out, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["genome1", "genome2", "chr1", "chr2", "startBp1", "endBp1", "startBp2", "endBp2",
                     "orient", "blkID", "refGenome", "refChr", "highlight"])
        wr.writerows(rows)
    for p in pairs:
        pr = [r for r in rows if r[0] == p["g1"] and r[1] == p["g2"]]
        bp = sum(r[5] - r[4] + 1 for r in pr if not r[12])
        lab_ = sum(1 for r in pr if r[11] and not r[12])
        info(f"{p['g1']} -> {p['g2']}: {sum(1 for r in pr if not r[12])} block(s), {bp / 1e6:.1f} Mb of "
             f"{p['g1']}, {lab_} traced to a {ref} chromosome; {sum(1 for r in pr if r[12])} interval row(s)")
    info(f"wrote {len(rows)} rows to {a.out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("genomes", help="row order and HAL names")
    g.add_argument("--indir", required=True, help="--outdir of 9a")
    g.add_argument("--species-tree"); g.add_argument("--tree-map")
    g.add_argument("--species-table"); g.add_argument("--only-listed", action="store_true")
    g.add_argument("--drop-species", default="")
    g.add_argument("--out", required=True)
    b = sub.add_parser("blk", help="halSynteny PSLs -> GENESPACE blk CSV")
    b.add_argument("--genomes", required=True, help="TSV from the genomes subcommand")
    b.add_argument("--psl-dir", required=True)
    b.add_argument("--region", help="9a region.tsv: the interval to highlight / zoom on")
    b.add_argument("--ref", help="reference code (default: from --region)")
    b.add_argument("--window-kb", type=float, help="zoom: keep the projected interval +/- this much")
    b.add_argument("--merge-kb", type=float, default=200.0,
                   help="projected pieces of the interval closer than this are one piece")
    b.add_argument("--min-block-kb", type=float, default=0.0)
    b.add_argument("--min-chr-mb", type=float, default=0.0, help="leave out contigs drawn shorter than this")
    b.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.cmd == "genomes" and a.only_listed and not a.species_table:
        sys.exit("[error] --only-listed needs --species-table")
    (cmd_genomes if a.cmd == "genomes" else cmd_blk)(a)


if __name__ == "__main__":
    main()
