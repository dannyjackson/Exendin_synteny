#!/usr/bin/env python3
"""
Rename tips in Newick trees and FASTA files with a two-column map (old<TAB>new), e.g.
5b's _map.tsv. Branch lengths and support values are untouched, so a tree already built
in step 4 is relabelled without re-running IQ-TREE.

Each input is written to --outdir under the same file name. Names not in the map are kept
and listed, so a missing name is seen, not silently left as the old one.

Usage
  python3 5c_rename_tips.py --map names_map.tsv --outdir renamed \
      ../reciprocal_blast/tree/exendin_aa.rooted.nwk ../reciprocal_blast/tree/exendin_aa.aln.fa
"""
import argparse, os, re, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", required=True, help="old<TAB>new")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("files", nargs="+", help="Newick (.nwk, .treefile, .tre) or FASTA files")
    a = ap.parse_args()
    m = {}
    for line in open(a.map):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[0] and not line.startswith("#"):
            m[f[0]] = f[1]
    os.makedirs(a.outdir, exist_ok=True)
    for path in a.files:
        txt = open(path).read()
        missing = set()
        if txt.lstrip().startswith(">"):
            def rep(mo):
                n = mo.group(1)
                if n not in m:
                    missing.add(n)
                return ">" + m.get(n, n) + mo.group(2)
            out = re.sub(r"^>(\S+)(.*)$", rep, txt, flags=re.M)
        else:
            def rep(mo):                       # a tip follows '(' or ','; node labels follow ')'
                n = mo.group(2).strip("'")
                if n not in m:
                    missing.add(n)
                return mo.group(1) + m.get(n, n)
            out = re.sub(r"([(,])\s*('[^']+'|[^(),:;]+)", rep, txt)
        dest = os.path.join(a.outdir, os.path.basename(path))
        open(dest, "w").write(out)
        print(f"[info] {path} -> {dest}" + (f"; not in the map, kept: {', '.join(sorted(missing))}"
                                           if missing else ""), file=sys.stderr)


if __name__ == "__main__":
    main()
