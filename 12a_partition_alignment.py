#!/usr/bin/env python3
"""
Split the curated exendin alignment into a SLOW partition (signal peptide + pro-region)
and a FAST partition (the mature peptide onward), so the topology can be built from the
slow part and the fast part mapped onto it (12_partitioned_tree.md).

Why. Irwin & Prentice (2011) found that on the Heloderma lineage the exendin toxin peptide
sequences, but not their pro or signal peptides, evolved very rapidly. A fast-evolving block
bolted onto a slow backbone is the standard way to drag a clade toward the root, which is
the likeliest explanation for the exendins sitting sister to all GCGL with weak support.
Partitioning separates the two signals: the slow block carries the phylogeny, the fast block
carries the adaptive story.

Boundary. The mature peptide start is found per sequence with the same rule as
5a_assign_nomenclature.py: the hormone-domain motif [HY][ASG][DEQ][AGP][ATVIGSM]F[TSI],
else the cleavage site (KR|RR|KK|R|K) nearest the start of a hormone-like window, else the
best-scoring window against the published panel. Each start is mapped to its alignment
column; the consensus column is the partition boundary. The spread across sequences is
reported -- if it is wide, the alignment is misaligned over the cleavage site and should be
fixed before trusting either partition.

Outputs (--outdir)
  boundaries.tsv       per sequence: how the mature peptide was found, its start as an
                       unaligned residue and as an alignment column, and the peptide
  slow.aa.fa, fast.aa.fa        the two protein partitions
  slow.codon.fa, fast.codon.fa  the same split of --codon-aln, if given
  partitions.nex                IQ-TREE/RAxML partition file for the codon alignment
  partitions_aa.nex             the same for the protein alignment

Usage
  python3 12a_partition_alignment.py --aln aln1/exendin.aa.fa --codon-aln aln1/exendin.codon.fa \
      --reference ${SCRIPTDIR}/4_reference_peptides.fa --outdir partitioned
"""
import argparse, importlib.util, os, statistics, sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("nomen", os.path.join(HERE, "5a_assign_nomenclature.py"))
NM = importlib.util.module_from_spec(spec)
spec.loader.exec_module(NM)


def read_aln(path):
    out, name = [], None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            out.append([name, []])
        elif out:
            out[-1][1].append(line.strip())
    recs = [(n, "".join(s).upper()) for n, s in out]
    if len({len(s) for _, s in recs}) != 1:
        sys.exit(f"[error] {path}: sequences differ in length, not an alignment")
    return recs


def col_of(aligned, residue_index):
    """1-based alignment column holding the given 0-based ungapped residue index."""
    k = -1
    for c, ch in enumerate(aligned):
        if ch not in "-.":
            k += 1
            if k == residue_index:
                return c + 1
    return None


def write_fa(path, recs, lo, hi):
    """columns lo..hi inclusive, 1-based."""
    with open(path, "w") as fh:
        for n, s in recs:
            fh.write(f">{n}\n{s[lo - 1:hi]}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--aln", required=True, help="curated protein alignment (aln1/exendin.aa.fa)")
    ap.add_argument("--codon-aln", help="the matching codon alignment (aln1/exendin.codon.fa)")
    ap.add_argument("--reference", required=True, help="4_reference_peptides.fa (published *_mature peptides)")
    ap.add_argument("--boundary", type=int, help="force the boundary: first alignment column of the FAST partition")
    ap.add_argument("--exclude", default="", help="comma-separated tips to ignore when choosing the consensus "
                                                  "(e.g. partial models, or GCG/GIP whose domains sit elsewhere)")
    ap.add_argument("--max-spread", type=int, default=6,
                    help="warn if the per-sequence boundary columns span more than this")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    aln = read_aln(a.aln)
    ncol = len(aln[0][1])
    refs = [(h.split("|")[0], s) for h, s in NM.read_fasta(a.reference) if h.split("|")[0].endswith("_mature")]
    if not refs:
        sys.exit(f"[error] no *_mature records in {a.reference}")
    skip = {x for x in a.exclude.split(",") if x}

    rows, cols = [], []
    for n, aligned in aln:
        full = aligned.replace("-", "").replace(".", "")
        pep, how = NM.mature(full, refs)
        start = full.find(pep) if pep else -1
        col = col_of(aligned, start) if start >= 0 else None
        rows.append([n, how, start + 1 if start >= 0 else "-", col or "-", len(pep), pep])
        if col and n not in skip:
            cols.append(col)
    if not cols:
        sys.exit("[error] no mature peptide start could be located in any sequence")

    counts = Counter(cols)
    consensus = counts.most_common(1)[0][0]
    spread = max(cols) - min(cols)
    boundary = a.boundary or consensus
    print(f"[info] mature-peptide start column: consensus {consensus} "
          f"(median {int(statistics.median(cols))}, range {min(cols)}-{max(cols)}, "
          f"{counts[consensus]}/{len(cols)} sequences agree)", file=sys.stderr)
    for c, k in sorted(counts.items()):
        if c != consensus:
            who = [r[0] for r in rows if r[3] == c and r[0] not in skip]
            print(f"[warn] {k} sequence(s) put it at column {c}: {', '.join(who[:6])}", file=sys.stderr)
    if spread > a.max_spread:
        print(f"[warn] boundary columns span {spread} columns (> --max-spread {a.max_spread}): the alignment "
              f"is probably ragged over the cleavage site. Fix it, or set --boundary by eye, before "
              f"trusting the partitions", file=sys.stderr)
    for r in rows:
        if r[1] != "motif":
            print(f"[warn] {r[0]}: mature peptide found by '{r[1]}', not the motif: check it", file=sys.stderr)
    if a.boundary:
        print(f"[info] --boundary {a.boundary} overrides the consensus", file=sys.stderr)
    if not 2 <= boundary <= ncol:
        sys.exit(f"[error] boundary {boundary} is outside the alignment (1-{ncol})")

    with open(os.path.join(a.outdir, "boundaries.tsv"), "w") as fh:
        fh.write("tip\thow_found\tmature_start_residue\tmature_start_column\tpeptide_len\tmature_peptide\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    slow = (1, boundary - 1)
    fast = (boundary, ncol)
    write_fa(os.path.join(a.outdir, "slow.aa.fa"), aln, *slow)
    write_fa(os.path.join(a.outdir, "fast.aa.fa"), aln, *fast)
    with open(os.path.join(a.outdir, "partitions_aa.nex"), "w") as fh:
        fh.write("#nexus\nbegin sets;\n"
                 f"  charset slow_signal_pro = {slow[0]}-{slow[1]};\n"
                 f"  charset fast_mature = {fast[0]}-{fast[1]};\n"
                 "end;\n")
    msg = (f"[info] protein alignment {ncol} columns -> slow {slow[0]}-{slow[1]} "
           f"({slow[1] - slow[0] + 1} aa), fast {fast[0]}-{fast[1]} ({fast[1] - fast[0] + 1} aa)")
    if a.codon_aln:
        cod = read_aln(a.codon_aln)
        ncod = len(cod[0][1])
        if ncod != ncol * 3:
            sys.exit(f"[error] {a.codon_aln} has {ncod} columns, not 3 x {ncol}: it does not match --aln")
        if {n for n, _ in cod} != {n for n, _ in aln}:
            sys.exit("[error] the codon and protein alignments hold different tips")
        cs, cf = (1, (boundary - 1) * 3), ((boundary - 1) * 3 + 1, ncod)
        write_fa(os.path.join(a.outdir, "slow.codon.fa"), cod, *cs)
        write_fa(os.path.join(a.outdir, "fast.codon.fa"), cod, *cf)
        with open(os.path.join(a.outdir, "partitions.nex"), "w") as fh:
            fh.write("#nexus\nbegin sets;\n"
                     f"  charset slow_signal_pro = {cs[0]}-{cs[1]};\n"
                     f"  charset fast_mature = {cf[0]}-{cf[1]};\n"
                     "end;\n")
        msg += f"; codon alignment -> slow {cs[0]}-{cs[1]}, fast {cf[0]}-{cf[1]}"
    print(msg, file=sys.stderr)
    print(f"[info] wrote {a.outdir}/boundaries.tsv and the partition files", file=sys.stderr)


if __name__ == "__main__":
    main()
