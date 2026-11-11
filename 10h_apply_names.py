#!/usr/bin/env python3
"""
Put the step-5 names on everything downstream: rewrite the tips of alignments and trees
from names/names.tsv (10f), so every later analysis and figure speaks one nomenclature
(<Sp>_GCGL, <Sp>_Exendin_anc, <Sp>_A/_B/_Z ..., and the unchanged GCG / GIP outgroups)
instead of the provisional TOGA2 labels (<Sp>_Ex2_Ex4b ...).

10f renames the tree it draws, but leaves the alignments and the step-4 helper files on the
old labels. Run this once after step 5; steps 6 and 12 then take the renamed files.

Also writes, in the new nomenclature:
  groups.tsv    tip -> class (GCGL, Exendin_anc, copy, GCG, GIP), for 4h_draw_tree.py --groups
  outgroup.txt  the GIP tips, comma-separated, for IQ-TREE -o
  map.tsv       old -> new, for anything renamed by hand later

A tip with no row in names.tsv keeps its name and is reported; nothing is silently dropped.

Usage
  python3 10h_apply_names.py --names names/names.tsv \
      --fasta aln1/exendin.aa.fa --fasta aln1/exendin.codon.fa \
      --tree tree/exendin_aa.treefile --outdir named
"""
import argparse, csv, os, re, sys


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--names", required=True, help="names/names.tsv from 10f")
    ap.add_argument("--fasta", action="append", default=[], help="alignment or FASTA to rename (repeatable)")
    ap.add_argument("--tree", action="append", default=[], help="Newick file to rename (repeatable)")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    rows = read_tsv(a.names)
    need = {"old_tip", "new_name", "class"}
    if not rows or not need <= set(rows[0]):
        sys.exit(f"[error] {a.names} needs the columns {', '.join(sorted(need))}: is it names/names.tsv from 10f?")
    mapping = {r["old_tip"]: r["new_name"] for r in rows}
    cls = {r["new_name"]: r["class"] for r in rows}
    dup = [n for n in set(mapping.values()) if list(mapping.values()).count(n) > 1]
    if dup:
        sys.exit(f"[error] {a.names} maps two tips to the same name: {', '.join(sorted(dup))}")
    print(f"[info] {len(mapping)} name(s) from {a.names}; "
          + ", ".join(f"{k} {sum(v == k for v in cls.values())}" for k in sorted(set(cls.values()))),
          file=sys.stderr)

    unmapped = set()

    def rename(tip):
        if tip in mapping:
            return mapping[tip]
        unmapped.add(tip)
        return tip

    for path in a.fasta:
        out = os.path.join(a.outdir, os.path.basename(path))
        n = 0
        with open(path) as fh, open(out, "w") as fo:
            for line in fh:
                if line.startswith(">"):
                    head = line[1:].rstrip("\n")
                    tip = head.split()[0]
                    rest = head[len(tip):]
                    fo.write(f">{rename(tip)}{rest}\n")
                    n += 1
                else:
                    fo.write(line)
        print(f"[info] {out}: {n} sequence(s)", file=sys.stderr)

    for path in a.tree:
        out = os.path.join(a.outdir, os.path.basename(path))
        txt = open(path).read()
        # tip labels sit after '(' or ',', before ':' or ')' or ','
        def sub(m):
            return m.group(1) + rename(m.group(2).strip("'"))
        new = re.sub(r"([(,])\s*('[^']*'|[^(),:;]+)", sub, txt)
        open(out, "w").write(new)
        print(f"[info] {out}", file=sys.stderr)

    with open(os.path.join(a.outdir, "groups.tsv"), "w") as fh:
        fh.write("tip\tgroup\n")
        for new, c in sorted(cls.items()):
            fh.write(f"{new}\t{c}\n")
    og = sorted(n for n, c in cls.items() if c == "GIP")
    with open(os.path.join(a.outdir, "outgroup.txt"), "w") as fh:
        fh.write(",".join(og) + "\n")
    with open(os.path.join(a.outdir, "map.tsv"), "w") as fh:
        fh.write("old_tip\tnew_name\n")
        for o, n in sorted(mapping.items()):
            fh.write(f"{o}\t{n}\n")
    if not og:
        print("[warn] no tip is classed GIP: outgroup.txt is empty, so IQ-TREE -o has nothing to root on",
              file=sys.stderr)
    if unmapped:
        print(f"[warn] {len(unmapped)} tip(s) have no row in {a.names} and keep their name: "
              f"{', '.join(sorted(unmapped)[:8])}", file=sys.stderr)
    print(f"[info] wrote {a.outdir}/groups.tsv, outgroup.txt ({len(og)} GIP tip(s)) and map.tsv",
          file=sys.stderr)


if __name__ == "__main__":
    main()
