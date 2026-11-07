#!/usr/bin/env python3
"""
Give every curated exendin-family gene (step 4) a standard name.

The names inherited from annotations (NCBI 'exendin-3-like', TOGA2 'exendin-2/3') are
neither monophyletic nor consistent between sources, so they are replaced. The working
hypothesis is one ancestral gene (GCGL) that duplicated several times; the rules follow
from it and make no orthology claim beyond it:

  1. one gene in a genome            -> <species>_GCGL
     (the ancestral single-copy state; every curated locus of the species counts,
     remnants and disabled copies included, so a genome with GCGL plus a pseudogene is
     not single-copy)
  2. several genes in a genome       -> <species>_A, _B, _C ... in order from SLC4A8 to AQP5
     Letters are positions, not orthology: gene A of one species need not be orthologous
     to gene A of another.
  3. Heloderma (--heloderma, default HelSus,HelCha)
                                     -> named by the published peptide its mature peptide
     matches (5a): EX1, EX2, EX3, EX4, or GCGL for the GCGL-type copy, when the identity
     is >= --close (80%). A copy that matches nothing published at that level keeps its
     positional letter. If two copies match the same peptide, the better match (e.g. the
     exact one) keeps the plain name and the other gets its positional letter appended
     (HelSus_EX4, HelSus_EX4_D); if they tie, both get letters.

Order from SLC4A8 to AQP5 (rule 2): on the contig carrying both boundary genes (4d
boundaries_used.tsv), loci are ordered by position measured from SLC4A8 towards AQP5,
whichever way round the assembly has them; loci outside the boundaries (HelSus has copies
beyond them) take their place on the same line. In a fragmented assembly: loci on the
SLC4A8 contig (nearest SLC4A8 first), then unplaced contigs, then loci on the AQP5
contig (furthest from AQP5 first), then 'elsewhere' loci. The 'order_basis' column says
which applied, so a letter that rests on an assembly guess can be seen.

Extra tips (e.g. the chicken outgroup added with 4g --outgroup-fasta) are not counted;
they are named <species>_GCGL if 5a calls them GCGL-type (>= --close) and keep their
name otherwise. If the same gene is also a curated locus (chicken's partial model), give
--drop to leave the duplicate out of the counts.

Inputs
  --curated     4d curated_loci.tsv
  --boundaries  4d boundaries_used.tsv
  --assignments 5a <prefix>_assignments.tsv, run on the tree's proteins (4g tree_aa.fa)
Outputs (--out-prefix)
  _names.tsv     old_name new_name species rule order_basis position placement n_in_genome
                 best_match identity verdict
  _map.tsv       old_name<TAB>new_name, for 5c_rename_tips.py
  _groups.tsv    new_name<TAB>colour class, for 4h --groups: the published peptide matched at
                 >= --close, else 'no close published match' (grey)
  _rules.tsv     new_name, code, detail, for 4h --sources (the code after each tip)
  _legend.txt    code=meaning list, for 4h --legend
"""
import argparse, csv, os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

# 5a best_match base -> short name (Heloderma rule) and colour class
SHORT = {"EXENDIN1": "EX1", "EXENDIN2": "EX2", "EXENDIN3": "EX3", "EXENDIN4": "EX4", "GCGL": "GCGL"}
CLASS = {"EXENDIN1": "exendin-1", "EXENDIN2": "exendin-2", "EXENDIN3": "exendin-3", "EXENDIN4": "exendin-4",
         "GCGL": "GCGL (glucagon-like)", "GCG": "proglucagon", "GLP1": "GLP-1", "GLP2": "GLP-2", "PACAP": "PACAP"}


def letters(i):
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--curated", required=True, help="4d curated_loci.tsv")
    ap.add_argument("--boundaries", required=True, help="4d boundaries_used.tsv")
    ap.add_argument("--assignments", required=True, help="5a _assignments.tsv")
    ap.add_argument("--heloderma", default="HelSus,HelCha", help="species named by published-peptide match")
    ap.add_argument("--close", type=float, default=80.0, help="identity for a published-peptide name (5a's 'close')")
    ap.add_argument("--drop", default="", help="comma-separated curated loci to leave out (duplicates of extra tips)")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    helo = set(x for x in a.heloderma.split(",") if x)
    drop = set(x for x in a.drop.split(",") if x)
    with open(a.curated) as fh:
        cur = [r for r in csv.DictReader(fh, delimiter="\t") if r["name"] not in drop]
    with open(a.boundaries) as fh:
        bnd = list(csv.DictReader(fh, delimiter="\t"))
    with open(a.assignments) as fh:
        asg = {r["gene"]: r for r in csv.DictReader(fh, delimiter="\t")}
    for n in drop:
        N.info(f"--drop {n}: left out of the counts and the names")

    def base(gene):
        r = asg.get(gene)
        return (r["best_match"].replace("_mature", ""), float(r["identity"])) if r else (None, 0.0)

    by_sp = defaultdict(list)
    for r in cur:
        by_sp[r["species"]].append(r)
    rows = []
    for sp, Ls in sorted(by_sp.items()):
        b = defaultdict(list)
        for x in bnd:
            if x["species"] == sp:
                b[x["boundary"]].append(x)
        s = b.get("SLC4A8", [None])[0]
        q = b.get("AQP5", [None])[0]
        mid = lambda x: (int(x["start"]) + int(x["end"])) / 2
        lmid = lambda r: (int(r["g_start"]) + int(r["g_end"])) / 2

        def key(r):
            if s and q and s["contig"] == q["contig"] and r["contig"] == s["contig"]:
                d = 1 if mid(q) > mid(s) else -1
                return (1, d * (lmid(r) - mid(s))), "interval, SLC4A8->AQP5"
            if s and r["contig"] == s["contig"]:
                return (1, abs(lmid(r) - mid(s))), "SLC4A8 contig, from SLC4A8"
            if q and r["contig"] == q["contig"]:
                return (3, -abs(lmid(r) - mid(q))), "AQP5 contig, towards AQP5"
            if r["placement"] == "unplaced" or not (s or q):
                return (2, r["contig"], lmid(r)), "unplaced contig: position assumed between the boundaries"
            return (4, r["contig"], lmid(r)), "elsewhere: not in the interval"
        Ls = sorted(Ls, key=lambda r: key(r)[0])
        n = len(Ls)
        names = []
        for i, r in enumerate(Ls):
            letter = letters(i)
            bm, ident = base(r["name"])
            if sp in helo:
                if bm in SHORT and ident >= a.close:
                    nm, rule = SHORT[bm], "H"
                else:
                    nm, rule = letter, "L"
            elif n == 1:
                nm, rule = "GCGL", "S"
            else:
                nm, rule = letter, "L"
            names.append([r, nm, rule, letter, key(r)[1], bm, ident])
        # a published name used twice: the copy matching it best keeps the plain name (if it is
        # strictly the best), the others get their positional letter appended
        dup = defaultdict(list)
        for x in names:
            if x[2] == "H":
                dup[x[1]].append(x)
        keep_plain = set()
        for nm, xs in dup.items():
            xs = sorted(xs, key=lambda x: -x[6])
            if len(xs) == 1 or xs[0][6] > xs[1][6]:
                keep_plain.add(id(xs[0]))
        for x in names:
            r, nm, rule, letter, basis, bm, ident = x
            if rule == "H" and len(dup[nm]) > 1 and id(x) not in keep_plain:
                nm = f"{nm}_{letter}"
            new = f"{sp}_{nm}"
            ar = asg.get(r["name"], {})
            rows.append(dict(old_name=r["name"], new_name=new, species=sp, rule=rule, order_basis=basis,
                             position=letter, placement=r["placement"], n_in_genome=n,
                             best_match=bm or "NA", identity=f"{ident:.1f}", verdict=ar.get("verdict", "not in 5a")))
        if sp in helo:
            N.info(f"{sp}: {n} gene(s): " + ", ".join(f"{x['new_name']}" for x in rows if x["species"] == sp))
        if any("assumed" in x[4] or "elsewhere" in x[4] for x in names):
            N.warn(f"{sp}: some letters rest on an unplaced or out-of-interval locus (see order_basis)")

    curated_names = {r["name"] for r in cur}
    for gene, r in asg.items():                       # extra tips (outgroup)
        if gene in curated_names or gene in drop:
            continue
        bm, ident = base(gene)
        sp = gene.split("_")[0]
        if gene.startswith("Published_"):               # 4g --published: keep the published name
            rows.append(dict(old_name=gene, new_name=gene, species="published", rule="P",
                             order_basis="published precursor, not counted", position="-", placement="-",
                             n_in_genome="-", best_match=bm or "NA", identity=f"{ident:.1f}",
                             verdict=r.get("verdict", "")))
            continue
        new = f"{sp}_GCGL" if bm == "GCGL" and ident >= a.close else gene
        rows.append(dict(old_name=gene, new_name=new, species=sp, rule="O", order_basis="extra tip, not counted",
                         position="-", placement="-", n_in_genome="-", best_match=bm or "NA",
                         identity=f"{ident:.1f}", verdict=r.get("verdict", "")))
        if new == gene:
            N.warn(f"{gene}: extra tip not GCGL-type at >= {a.close:g}% ({bm} {ident:.0f}%); name kept")
    seen = defaultdict(list)
    for x in rows:
        seen[x["new_name"]].append(x["old_name"])
    for k, v in seen.items():
        if len(v) > 1:
            N.warn(f"new name {k} given to {', '.join(v)}: rename by hand in _map.tsv")

    cols = ["old_name", "new_name", "species", "rule", "order_basis", "position", "placement", "n_in_genome",
            "best_match", "identity", "verdict"]
    with open(a.out_prefix + "_names.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for x in rows:
            fh.write("\t".join(str(x[c]) for c in cols) + "\n")
    with open(a.out_prefix + "_map.tsv", "w") as fh:
        for x in rows:
            fh.write(f"{x['old_name']}\t{x['new_name']}\n")
    with open(a.out_prefix + "_groups.tsv", "w") as fh:
        fh.write("tip\tgroup\n")
        for x in rows:
            if x["best_match"] == "NA":
                cls = "not compared"
            elif float(x["identity"]) < a.close:
                cls = "no close published match"     # most squamate copies: ~50% to anything published
            else:
                cls = CLASS.get(x["best_match"], x["best_match"].lower())
            fh.write(f"{x['new_name']}\t{cls}\n")
    with open(a.out_prefix + "_rules.tsv", "w") as fh:
        fh.write("tip\tcode\tsource\tspecies\tprotein\tdescription\n")
        for x in rows:
            ident = f"{x['best_match']} {x['identity']}%"
            code = {"H": f"H {x['identity'].split('.')[0]}%", "L": "L", "S": "S", "O": "O", "P": "P"}[x["rule"]]
            fh.write(f"{x['new_name']}\t{code}\t{x['old_name']}\t{x['species']}\t{x['order_basis']}\t"
                     f"closest published: {ident}\n")
    with open(a.out_prefix + "_legend.txt", "w") as fh:
        fh.write(f"S=only gene in its genome: GCGL\nL=position from SLC4A8 to AQP5 (not orthology)\n"
                 f"H=Heloderma: published peptide matched, % identity (>= {a.close:g}%)\nO=outgroup\n"
                 f"P=published precursor (UniProt)\n")
    k = defaultdict(int)
    for x in rows:
        k[x["rule"]] += 1
    N.info(f"{len(rows)} genes named: " + ", ".join(f"{c} {k[c]}" for c in "SLHOP" if k[c])
           + f"; written to {a.out_prefix}_*.tsv")


if __name__ == "__main__":
    main()
