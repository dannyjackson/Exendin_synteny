#!/usr/bin/env python3
"""
Step 14: WormSynteny-style tables (seqs, genes, links) for gggenomes, from the two halLiftover
waves of 14_wormsynteny.md.

WormSynteny (Zheng lab; wormsynteny.org) draws one region of a Cactus alignment as follows:
  1. halLiftover the query region into every species (first wave),
  2. stitch fragments closer than a gap value; drop fragments shorter than a filtering
     percentage of the query (bedtools merge + awk in 14_wormsynteny.md),
  3. keep the longest fragment per species as its row; the others are "(filtered)" extras,
  4. halLiftover each species' fragments into the next species of the row order (second
     wave): these alignments are the links,
  5. gggenomes(seqs, genes, links), genes coloured by orthogroup.
Its package is hard-wired to its own Caenorhabditis HAL, species and annotation, so this
script applies the same steps to this project's HAL: the query is the HelSus SLC4A8-AQP5
interval (9a region.tsv), the genes are each genome's own TOGA2 annotation (9a
target_genes.tsv / ref_genes.tsv), and the "orthogroup" of a gene is the HelSus gene TOGA2
projected it from: the step-5 name for exendin-family genes (anc, A, B, Z ...), else the symbol.

Outputs (<prefix>_seqs.tsv, _genes.tsv, _links.tsv) are in gggenomes' column names, with
every coordinate relative to its fragment (1 = fragment start). Standard library only.

Usage
  python3 14a_wormsynteny_tables.py --genomes ../genespace_cactus/genomes.tsv \
      --lift ../synteny_cactus/hal_lift --frag-dir frags --pair-dir pairs \
      --gap 5000 --out wormsynteny 2> wormsynteny_tables.log
"""
import argparse, csv, os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
SHORT = {"exendin_anc": "anc", "gcgl": "GCGL"}


def info(m):
    print(f"[info] {m}", file=sys.stderr)


def warn(m):
    print(f"[warn] {m}", file=sys.stderr)


def read_tsv(path):
    if not os.path.exists(path):
        sys.exit(f"[error] no such file: {path}")
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def canon(c):
    """NC_086153.1 = NC_086153__1 = NC_086153_1 (TOGA2 and Cactus spell the version differently)."""
    return re.sub(r"[._]+", "_", (c or "").strip().lower())


def read_bed(path):
    out = []
    if os.path.exists(path):
        for line in open(path):
            t = line.rstrip("\n").split("\t")
            if len(t) >= 3 and not line.startswith(("#", "track")):
                out.append((t[0], int(t[1]) + 1, int(t[2])))
    return out


def psl_blocks(path):
    """halLiftover --outPSL: q = source genome, t = destination. -> forward-strand blocks
    (qchr, qs, qe, tchr, ts, te, strand), 1-based inclusive."""
    out = []
    if not os.path.exists(path):
        return out
    for line in open(path):
        t = line.rstrip("\n").split("\t")
        if len(t) < 21 or not t[0].isdigit():
            continue
        qst, tst = t[8][0], t[8][-1]
        qsize, tsize = int(t[10]), int(t[14])
        sizes = [int(x) for x in t[18].rstrip(",").split(",")]
        qss = [int(x) for x in t[19].rstrip(",").split(",")]
        tss = [int(x) for x in t[20].rstrip(",").split(",")]
        for n, qs, ts in zip(sizes, qss, tss):
            if qst == "-":
                qs = qsize - qs - n
            if tst == "-":
                ts = tsize - ts - n
            out.append((t[9], qs + 1, qs + n, t[13], ts + 1, ts + n, "+" if qst == tst else "-"))
    return out


def chain(blocks, gap):
    """blocks on one fragment pair -> links: consecutive blocks (in query order) on the same
    strand, with both gaps <= gap, become one link (as halSynteny / WormSynteny stitching)."""
    out = []
    for b in sorted(blocks, key=lambda b: (b[6], b[1])):
        qs, qe, ts, te, st = b[1], b[2], b[4], b[5], b[6]
        c = out[-1] if out else None
        if c and c["st"] == st and 0 <= qs - c["qe"] - 1 <= gap and (
                (st == "+" and -gap <= ts - c["te"] - 1 <= gap) or (st == "-" and -gap <= c["ts"] - te - 1 <= gap)):
            c["qe"] = max(c["qe"], qe); c["ts"] = min(c["ts"], ts); c["te"] = max(c["te"], te)
        else:
            out.append(dict(qs=qs, qe=qe, ts=ts, te=te, st=st))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--genomes", required=True, help="genomes.tsv of 13a (row order, HAL names)")
    ap.add_argument("--lift", required=True, help="9a --outdir: region.tsv, ref_genes.tsv, target_genes.tsv")
    ap.add_argument("--frag-dir", required=True, help="first wave, stitched and filtered: <hal>.bed")
    ap.add_argument("--pair-dir", required=True, help="second wave: <halA>__<halB>.psl")
    ap.add_argument("--species-table", help="species.tsv, for the row labels")
    ap.add_argument("--gap", type=int, default=5000, help="links closer than this are chained (bp)")
    ap.add_argument("--extras", action="store_true",
                    help="also draw each genome's other fragments, as '(filtered)' rows at the bottom")
    ap.add_argument("--min-link-bp", type=int, default=200, help="links shorter than this are not drawn")
    ap.add_argument("--out", required=True, help="output prefix")
    a = ap.parse_args()

    g = read_tsv(a.genomes)
    order = [r["code"] for r in g]
    halof = {r["code"]: r["hal_genome"] for r in g}
    reg = read_tsv(os.path.join(a.lift, "region.tsv"))[0]
    lab = {}
    if a.species_table:
        import lib_exendin_names as N
        lab = {r["code"]: r["species"].replace("_", " ") for r in N.read_species_table(a.species_table)}

    # ---- seqs: the longest fragment of each genome is its row ----
    main, extras = {}, defaultdict(list)
    for c in order:
        fr = read_bed(os.path.join(a.frag_dir, f"{halof[c]}.bed"))
        if not fr:
            warn(f"{c}: no fragment passed the gap / filtering step; row left out")
            continue
        fr.sort(key=lambda f: f[2] - f[1], reverse=True)
        main[c], extras[c] = fr[0], fr[1:]
        info(f"{c}: row = {fr[0][0]}:{fr[0][1]}-{fr[0][2]} ({(fr[0][2] - fr[0][1] + 1) / 1000:.0f} kb)"
             + (f"; {len(fr) - 1} other fragment(s)" + ("" if a.extras else ", not drawn (--extras)") if len(fr) > 1 else ""))
    rows = [c for c in order if c in main]
    sid = lambda c, f: f"{c} {f[0]}:{f[1]}-{f[2]}"
    seqs = []
    for c in rows:
        f = main[c]
        seqs.append([sid(c, f), f[2] - f[1] + 1, lab.get(c, c), c, f[0], f[1], f[2], "main"])
    if a.extras:
        for c in rows:
            for f in extras[c]:
                seqs.append([sid(c, f), f[2] - f[1] + 1, f"{lab.get(c, c)} (filtered)", c, f[0], f[1], f[2], "filtered"])

    # ---- genes: each genome's own TOGA2 annotation on its fragments ----
    import lib_exendin_names as N
    tg = defaultdict(list)
    for r in read_tsv(os.path.join(a.lift, "target_genes.tsv")):
        tg[r["code"]].append(r)
    refc = order[0] if reg["ref"] not in halof.values() else next(c for c in order if halof[c] == reg["ref"])
    for r in read_tsv(os.path.join(a.lift, "ref_genes.tsv")):
        tg[refc].append(dict(r, code=refc, contig=reg["chrom"]))
    keys, genes = {}, []
    for s in seqs:
        seq_id, c, chrom, fs, fe = s[0], s[3], s[4], s[5], s[6]
        for r in tg.get(c, []):
            if canon(r["contig"]) != canon(chrom):
                continue
            gs, ge = int(r["start"]), int(r["end"])
            if ge < fs or gs > fe:
                continue
            name = r.get("name") or ""
            stem = name.partition("_")[2] if name.startswith(c + "_") else name
            stem = SHORT.get(stem.lower(), re.sub(r"^Exendin_", "", stem))
            og = stem if name else (r.get("symbol") or "").split("/")[0] or "unnamed"
            cov = "FULL" if gs >= fs and ge <= fe else "PARTIAL"
            genes.append([seq_id, max(gs, fs) - fs + 1, min(ge, fe) - fs + 1, r.get("strand") or "+",
                          r.get("gene_id", ""), og, og if name else "", int(bool(name)),
                          r.get("status", ""), cov, f"{r['contig']}:{gs}-{ge}"])
            keys[og] = keys.get(og) or bool(name)
    named = sorted([k for k, v in keys.items() if v])
    cmap = N.colour_map(named) if named else {}
    for gr in genes:
        gr.append(cmap.get(gr[5], ""))
    nown = defaultdict(int)
    for gr in genes:
        nown[gr[0].split(" ")[0]] += 1
    for c in rows:
        if not nown.get(c):
            warn(f"{c}: no TOGA2 gene on its fragment(s)")

    # ---- links: second wave, chained, between neighbouring rows' main fragments ----
    links = []
    for up, dn in zip(rows, rows[1:]):
        fa, fb = main[up], main[dn]
        path = os.path.join(a.pair_dir, f"{halof[up]}__{halof[dn]}.psl")
        if not os.path.exists(path):
            warn(f"{up} -> {dn}: no {path}; no links between these rows")
            continue
        bl = [b for b in psl_blocks(path)
              if b[0] == fa[0] and b[3] == fb[0] and b[2] >= fa[1] and b[1] <= fa[2] and b[5] >= fb[1] and b[4] <= fb[2]]
        ch = [l for l in chain(bl, a.gap) if l["qe"] - l["qs"] + 1 >= a.min_link_bp]
        for l in ch:
            links.append([sid(up, fa), max(l["qs"], fa[1]) - fa[1] + 1, min(l["qe"], fa[2]) - fa[1] + 1,
                          sid(dn, fb), max(l["ts"], fb[1]) - fb[1] + 1, min(l["te"], fb[2]) - fb[1] + 1, l["st"]])
        info(f"{up} -> {dn}: {len(bl)} aligned block(s) between the two rows -> {len(ch)} link(s), "
             f"{sum(l['qe'] - l['qs'] + 1 for l in ch) / 1000:.0f} kb of {up}")

    def w(name, head, rws):
        with open(f"{a.out}_{name}.tsv", "w") as fh:
            fh.write("\t".join(head) + "\n")
            for r in rws:
                fh.write("\t".join(map(str, r)) + "\n")
    w("seqs", ["seq_id", "length", "bin_id", "code", "contig", "contig_start", "contig_end", "type"], seqs)
    w("genes", ["seq_id", "start", "end", "strand", "gene_id", "og", "label", "named", "status", "coverage",
                "toga2_position", "colour"], genes)
    w("links", ["seq_id", "start", "end", "seq_id2", "start2", "end2", "strand"], links)
    info(f"wrote {a.out}_seqs.tsv ({len(seqs)}), _genes.tsv ({len(genes)}), _links.tsv ({len(links)})")


if __name__ == "__main__":
    main()
