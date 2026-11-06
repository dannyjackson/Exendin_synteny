#!/usr/bin/env python3
"""
List the genes that the HelSus reference annotation calls an exendin, and check
that each one is a reference gene in the TOGA2 run.

This defines the gene set for the TOGA2-only exendin tree (4_toga2_exendin_tree.md).
Nothing is added by hand and nothing comes from BLAST: a gene is in the set if and
only if its gene or transcript record in the reference GFF matches --pattern
(case-insensitive) in its attribute column.

Read-only on the GFF and the TOGA2 directory. Writes one table (--out, or stdout).

Usage
  python3 4a_list_reference_exendins.py \
      --gff ${HELSUS_EVIANN_GFF} \
      --run ${TOGADIR}/target_HelSus \
      --out reference_exendins.tsv

Output columns
  gene_id         reference gene ID (as used by TOGA2, e.g. LOC_00012325)
  label           short tip label for the tree (EX2, EX4, ...; EXN if no number in the
                  annotation). Edit this column in your copy of the table to rename tips;
                  4b reads it from there.
  chrom, start, end, strand
  transcripts     reference transcripts of the gene, from the GFF
  matched_text    the attribute(s) that matched --pattern, so you can see why it is in
  toga2_queries   queries whose orthology_classification.tsv lists the gene
  toga2_orthologs queries in which TOGA2 assigned it at least one ortholog
"""
import argparse, gzip, os, re, sys
from collections import defaultdict
from urllib.parse import unquote


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def attrs(col9):
    out = {}
    for kv in col9.strip().strip(";").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = unquote(v.strip())
    return out


def suggest_label(text):
    nums = sorted(set(re.findall(r"exendin[\s_-]*(\d+)", text, re.I)))
    return "EX" + "_".join(nums) if nums else "EXN"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True, help="the HelSus annotation TOGA2 projected (EviAnn GFF3)")
    ap.add_argument("--pattern", default="exendin", help="regex matched against gene/transcript attributes")
    ap.add_argument("--run", help="${TOGADIR}/target_HelSus, to confirm each gene is a TOGA2 reference gene")
    ap.add_argument("--out", help="output TSV (default stdout)")
    a = ap.parse_args()

    pat = re.compile(a.pattern, re.I)
    genes, tx_of, parent_of, hits = {}, defaultdict(list), {}, defaultdict(list)
    with opener(a.gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 9:
                continue
            ftype, at = t[2], attrs(t[8])
            fid = at.get("ID", "")
            if ftype == "gene":
                genes[fid] = (t[0], int(t[3]), int(t[4]), t[6])
            elif ftype in ("mRNA", "transcript"):
                par = at.get("Parent", "").split(",")[0]
                parent_of[fid] = par
                tx_of[par].append(fid)
            else:
                continue
            for k, v in at.items():
                if k not in ("ID", "Parent") and pat.search(v):
                    hits[fid].append(f"{ftype}:{k}={v}")

    sel = defaultdict(list)
    for fid, txt in hits.items():
        g = fid if fid in genes else parent_of.get(fid, fid)
        sel[g].extend(txt)
    if not sel:
        sys.exit(f"[error] no gene or transcript attribute in {a.gff} matches /{a.pattern}/")

    seen, present = defaultdict(set), defaultdict(set)
    if a.run:
        for qd in sorted(d for d in os.listdir(a.run) if d.startswith("query_")):
            p = os.path.join(a.run, qd, "orthology_classification.tsv")
            if not os.path.exists(p):
                print(f"[warn] missing {p}", file=sys.stderr)
                continue
            with opener(p) as fh:
                next(fh, None)
                for line in fh:
                    t = line.rstrip("\n").split("\t")
                    if len(t) >= 5 and t[0] in sel:
                        seen[t[0]].add(qd[6:])
                        if t[3] != "None":
                            present[t[0]].add(qd[6:])

    out = open(a.out, "w") if a.out else sys.stdout
    cols = ["gene_id", "label", "chrom", "start", "end", "strand", "transcripts",
            "matched_text", "toga2_queries", "toga2_orthologs"]
    print("\t".join(cols), file=out)
    used = defaultdict(int)
    for g in sorted(sel, key=lambda g: (genes.get(g, ("~", 0))[0], genes.get(g, ("~", 0))[1])):
        chrom, s, e, st = genes.get(g, ("NA", -1, -1, "."))
        text = " | ".join(dict.fromkeys(sel[g]))
        lab = suggest_label(text)
        used[lab] += 1
        if used[lab] > 1:
            lab = f"{lab}{'abcdefghijklmnopqrstuvwxyz'[used[lab] - 2]}"
        print("\t".join(map(str, [g, lab, chrom, s, e, st, ",".join(tx_of.get(g, [])) or "NA", text,
                                  len(seen[g]) if a.run else "NA",
                                  len(present[g]) if a.run else "NA"])), file=out)
        if a.run and not seen[g]:
            print(f"[warn] {g} matches /{a.pattern}/ but is in no orthology_classification.tsv: "
                  f"not a reference gene in this TOGA2 run", file=sys.stderr)
    if a.out:
        out.close()
    print(f"[info] {len(sel)} reference genes match /{a.pattern}/", file=sys.stderr)


if __name__ == "__main__":
    main()
