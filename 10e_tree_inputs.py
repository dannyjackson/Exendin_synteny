#!/usr/bin/env python3
"""
Tree inputs for the curated exendin alignment (10_toga2_exendin_tree.md, step 4): tip
colour groups, the outgroup, and constraint trees for testing where the exendins come from.

Tips are classified by name (after curation):
  *_GIP                     GIP          the outgroup (root)
  *_GCG                     proglucagon
  *_GCGL                    GCGL         GalGal_GCGL = chicken; the rest = squamate GCGL
  anything else             exendin      colour group from the name: Ex2 -> exendin-2,
                                         Ex2_Ex4 -> exendin-2/4, Ex3 -> exendin-3 ...
A tip whose name does not follow <Code>_<label> is reported and left out of the constraints.

Constraint trees (IQ-TREE -g; tips not named are free), one per hypothesis for the origin
of the exendins (all exendin tips together, 'EX'):
  H1_squamate_GCGL   EX + squamate GCGL form a clade that excludes chicken GCGL:
                     the exendins arose by duplication of GCGL within squamates
  H2_sister_GCGL     EX is sister to all GCGL (chicken included), together excluding GCG:
                     the duplication is older than the bird-squamate split
  H3_GCG             EX + GCG form a clade that excludes GCGL: the exendins derive from
                     proglucagon, not GCGL

Writes to --outdir: groups.tsv, outgroup.txt (comma-separated GIP tips), H1/H2/H3 .constraint.nwk,
tip_classes.tsv.

Usage
  python3 10e_tree_inputs.py --aln aln1/exendin.aa.fa --outdir tree
"""
import argparse, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402


def tips_of(path):
    return [l[1:].split()[0] for l in open(path) if l.startswith(">")]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--aln", required=True, help="the alignment the tree is built from")
    ap.add_argument("--chicken", default="GalGal", help="species code of the non-squamate GCGL")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    cls = {}
    for t in tips_of(a.aln):
        if "_" not in t:
            N.warn(f"{t}: not <Code>_<label>; left out of groups and constraints")
            continue
        code, label = t.split("_", 1)
        base = re.sub(r"(?<=[0-9A-Z])[b-z]$|p$", "", label)
        if base == "GIP":
            cls[t] = ("GIP", "GIP")
        elif base == "GCG":
            cls[t] = ("GCG", "proglucagon")
        elif base == "GCGL":
            cls[t] = ("GCGL_chicken" if code == a.chicken else "GCGL_squamate", "GCGL (glucagon-like)")
        else:
            g = N.name_group(base.upper())
            cls[t] = ("EX", g if g.startswith("exendin") else "exendin (unnumbered)")
    by = {k: [t for t, (c, _) in cls.items() if c == k] for k in ("GIP", "GCG", "GCGL_chicken", "GCGL_squamate", "EX")}
    for k, v in by.items():
        N.info(f"{k}: {len(v)} tip(s)" + (f": {', '.join(v)}" if len(v) <= 8 else ""))
    missing = [k for k in ("GIP", "GCG", "GCGL_chicken", "GCGL_squamate", "EX") if not by[k]]
    if missing:
        N.warn(f"no tips for {', '.join(missing)}: the constraint trees that need them are not written")

    with open(os.path.join(a.outdir, "groups.tsv"), "w") as fh:
        fh.write("tip\tgroup\n")
        for t, (_, g) in cls.items():
            fh.write(f"{t}\t{g}\n")
    with open(os.path.join(a.outdir, "tip_classes.tsv"), "w") as fh:
        fh.write("tip\tclass\tgroup\n")
        for t, (c, g) in cls.items():
            fh.write(f"{t}\t{c}\t{g}\n")
    with open(os.path.join(a.outdir, "outgroup.txt"), "w") as fh:
        fh.write(",".join(by["GIP"]) + "\n")

    j = ",".join
    EX, SQ, CH, GCG, GIP = by["EX"], by["GCGL_squamate"], by["GCGL_chicken"], by["GCG"], by["GIP"]
    out = lambda *xs: "".join("," + j(x) for x in xs if x)     # tips held outside the constrained clade
    cons = {}
    if EX and SQ and CH:
        cons["H1_squamate_GCGL"] = f"(({j(EX + SQ)}){out(CH, GCG, GIP)});"
    if EX and SQ + CH and GCG:
        cons["H2_sister_GCGL"] = f"((({j(EX)}),({j(SQ + CH)})){out(GCG, GIP)});"
    if EX and GCG and SQ + CH:
        cons["H3_GCG"] = f"(({j(EX + GCG)}){out(SQ + CH, GIP)});"
    for name, nwk in cons.items():
        with open(os.path.join(a.outdir, f"{name}.constraint.nwk"), "w") as fh:
            fh.write(nwk + "\n")
    N.info(f"wrote groups.tsv, outgroup.txt ({len(by['GIP'])} GIP tips) and {len(cons)} constraint tree(s): "
           + ", ".join(cons))


if __name__ == "__main__":
    main()
