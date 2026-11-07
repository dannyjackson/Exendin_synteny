#!/usr/bin/env python3
"""
Tie each curated gene to the published exendin nomenclature.

Copy of 2d_assign_nomenclature.py from Exendin_evolution with ONE fix: the old output line
wrote a literal '...' where the differences column belongs, so every row had one column
fewer than the header (differences and second_match ran together). Nothing else changed.

The names exendin-1 to -4 were given to peptides purified from venom, before any
genome existed, so they are defined by mature-peptide sequence and not by locus.
This script therefore compares the MATURE PEPTIDE of each curated protein against
the published peptides (2_reference_peptides.fa) by global alignment, and says for
each gene whether it encodes a published peptide exactly, differs from one by a few
residues, or matches nothing published.

The mature peptide is taken from the first residue of the hormone domain: the match
to [HY][ASG][DEQ][AGP][ATVIGSM]F[TSI], else the dibasic/monobasic cleavage site
(KR|RR|KK|R|K) closest to the start of a hormone-like window, else the best-scoring
window against the panel. The method used is reported for every gene; check any row
that is not 'motif' by eye.

Outputs (--out-prefix)
  _assignments.tsv   one row per gene: mature peptide, best published match,
                     identity, differences (ref->query at mature-peptide positions),
                     and a verdict
  _identity.tsv      full identity matrix, genes x published peptides

Usage
  python3 2d_assign_nomenclature.py --query curated/curated_aa.fa \
      --reference 2_reference_peptides.fa --out-prefix nomenclature
"""
import argparse, re, sys

MOTIF = re.compile(r"[HY][ASG][DEQ][AGP][ATVIGSM]F[TSI]")
CLEAV = re.compile(r"(KR|RR|KK|R|K)")
DIBASIC = re.compile(r"(KR|RR|KK)")

def read_fasta(path):
    out, name, buf = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:], []
        elif line:
            buf.append(line)
    if name:
        out.append((name, "".join(buf)))
    return out


def nw(a, b, match=1, mis=-1, gap=-2):
    """Global alignment; returns the two aligned strings."""
    n, m = len(a), len(b)
    S = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        S[i][0] = i * gap
    for j in range(1, m + 1):
        S[0][j] = j * gap
    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            S[i][j] = max(S[i - 1][j - 1] + (match if ai == b[j - 1] else mis),
                          S[i - 1][j] + gap, S[i][j - 1] + gap)
    i, j, A, B = n, m, [], []
    while i or j:
        if i and j and S[i][j] == S[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mis):
            A.append(a[i - 1]); B.append(b[j - 1]); i -= 1; j -= 1
        elif i and S[i][j] == S[i - 1][j] + gap:
            A.append(a[i - 1]); B.append("-"); i -= 1
        else:
            A.append("-"); B.append(b[j - 1]); j -= 1
    return "".join(reversed(A)), "".join(reversed(B))


def compare(query, ref):
    """Identity over aligned columns, plus differences numbered on the reference."""
    qa, ra = nw(query, ref)
    same = sum(1 for x, y in zip(qa, ra) if x == y and x != "-")
    cols = sum(1 for x, y in zip(qa, ra) if not (x == "-" and y == "-"))
    diffs, pos = [], 0
    for x, y in zip(qa, ra):
        if y != "-":
            pos += 1
        if x != y:
            diffs.append(f"{y if y != '-' else '-'}{pos}{x if x != '-' else '-'}")
    return 100.0 * same / cols if cols else 0.0, diffs

def trim(pep, minlen):
    """Cut the peptide at the first downstream dibasic processing site, if any."""
    m = DIBASIC.search(pep, minlen)
    return pep[:m.start()] if m else pep

def mature(seq, panel, minlen=24):
    """(peptide, how_it_was_found)"""
    m = MOTIF.search(seq)
    if m and len(seq) - m.start() >= minlen:
        return trim(seq[m.start():], minlen), "motif"
    best = None
    for c in CLEAV.finditer(seq):
        s = c.end()
        if len(seq) - s < minlen:
            continue
        sc = max(compare(seq[s:s + 40], p)[0] for _, p in panel)
        if seq[s] in "HY":
            sc += 10          # every peptide in this family starts with His or Tyr
        if best is None or sc > best[0]:
            best = (sc, s)
    win = None
    for s in range(0, max(1, len(seq) - minlen + 1)):
        sc = max(compare(seq[s:s + 40], p)[0] for _, p in panel)
        if win is None or sc > win[0]:
            win = (sc, s)
    # a real precursor is cut at a basic site, so keep a cleavage start unless a
    # free window beats it clearly
    if best and (win is None or best[0] >= win[0] - 5):
        return trim(seq[best[1]:], minlen), "cleavage_site"
    return (trim(seq[win[1]:], minlen), "best_window") if win else (seq, "whole_sequence")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", required=True, help="protein FASTA (step 5: tree/tree_aa.fa from 4g)")
    ap.add_argument("--reference", required=True, help="2_reference_peptides.fa")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--exact", type=float, default=100.0, help="identity counted as an exact match")
    ap.add_argument("--close", type=float, default=80.0, help="identity counted as the same peptide class")
    a = ap.parse_args()

    refs = [(h.split("|")[0], s) for h, s in read_fasta(a.reference) if h.split("|")[0].endswith("_mature")]
    meta = {h.split("|")[0]: h for h, _ in read_fasta(a.reference)}
    if not refs:
        sys.exit("no *_mature records in the reference file")
    queries = read_fasta(a.query)

    with open(a.out_prefix + "_assignments.tsv", "w") as out, \
         open(a.out_prefix + "_identity.tsv", "w") as mat:
        out.write("gene\tmature_peptide\tpeptide_len\tn_hormone_domains\thow_found\tbest_match\tidentity\t"
                  "differences\tsecond_match\tsecond_identity\tverdict\tbest_match_provenance\n")
        mat.write("gene\t" + "\t".join(n for n, _ in refs) + "\n")
        for h, seq in queries:
            gene = h.split()[0]
            full = seq.upper()
            pep, how = mature(full, refs)
            ndom, last = 0, -99
            for m in MOTIF.finditer(full):       # non-overlapping hormone domains
                if m.start() >= last + 20:
                    ndom += 1; last = m.start()
            def best_of(r):
                # an unprocessed C-terminal extension should not count as mismatch
                return max(compare(pep, r), compare(pep[:len(r)], r), key=lambda x: x[0])
            scored = sorted(((best_of(r)[0], best_of(r)[1], n) for n, r in refs), reverse=True)
            (id1, d1, n1), (id2, _, n2) = scored[0], scored[1]
            base = n1.replace("_mature", "")
            if id1 >= a.exact:
                verdict = f"encodes published {base} exactly"
            elif id1 >= a.close:
                verdict = f"{base}-type, {len(d1)} difference(s) from the published peptide"
            else:
                verdict = f"no published match (closest {base}, {id1:.0f}%)"
            out.write(f"{gene}\t{pep}\t{len(pep)}\t{ndom}\t{how}\t{n1}\t{id1:.1f}\t{','.join(d1) or '-'}\t"
                      f"{n2}\t{id2:.1f}\t{verdict}\t{meta[n1]}\n")
            mat.write(gene + "\t" + "\t".join(f"{compare(pep, r)[0]:.1f}" for _, r in refs) + "\n")
    print(f"[info] {len(queries)} genes compared with {len(refs)} published peptides", file=sys.stderr)


if __name__ == "__main__":
    main()
