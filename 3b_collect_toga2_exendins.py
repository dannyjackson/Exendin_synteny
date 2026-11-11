#!/usr/bin/env python3
"""
One sequence per TOGA2-annotated exendin gene in every query, for the TOGA2-only
exendin tree (4_toga2_exendin_tree.md).

The gene set is the table from 4a_list_reference_exendins.py: the genes the HelSus
reference annotation calls an exendin. Everything else comes from TOGA2 itself:

  - a tip is a TOGA2 query gene (q_gene in orthology_classification.tsv) that TOGA2
    assigns as an ortholog of one or more of those reference genes. HelSus is
    included through its self-projection (query_HelSus);
  - its sequence is TOGA2's own projected CDS (nucleotide.fa.gz), translated here.
    For a q_gene with several projections (isoforms, chains), the one with the best
    TOGA2 status is kept, then the lowest mRNA number, then the longest CDS;
  - the only filters are the TOGA2 status (--status) and a minimum CDS (--min-cds).

There is no curation file, no rescued model and no BLAST. Paralogous projections
(projections TOGA2 did not classify as orthologs) are left out unless
--include-paralogs is given; they are then merged into loci by overlap, loci that
overlap an ortholog are dropped, and tips get a 'p' suffix.

Read-only on the TOGA2 directory. Writes only to --outdir.

Usage
  python3 4b_collect_toga2_exendins.py \
      --run ${TOGADIR}/target_HelSus \
      --gene-table reference_exendins.tsv \
      --outdir toga2_exendins

Outputs (--outdir)
  projections.tsv   every projection of the reference genes, with status and class
  tips.tsv          one row per q_gene (and paralog locus): tip name, orthology,
                    chosen projection, status, coordinates, QC, included yes/no
  exendin_cds.fa    CDS of included tips (final stop removed)
  exendin_aa.fa     their translations
  groups.tsv        tip -> reference gene label(s), for 2c_draw_tree.py --groups

Tip names: <Query>_<label>, label(s) of the reference gene(s) the tip is orthologous
to (EX2_EX4 for a many2many gene). A letter suffix (b, c ...) separates several genes
with the same orthology in one species; tips.tsv maps each tip to its TOGA2 q_gene.
"""
import argparse, csv, gzip, os, re, sys
from collections import defaultdict

RANK = {"FI": 0, "I": 1, "PI": 2, "UL": 3, "M": 4, "L": 5, "PG": 6, "N": 7}
_B, _AA = "TCAG", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: _AA[16 * i + 4 * j + k]
        for i, a in enumerate(_B) for j, b in enumerate(_B) for k, c in enumerate(_B)}


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def translate(nt):
    return "".join(CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def gene_of(proj_or_tx):
    return proj_or_tx.split("-mRNA")[0]


def isoform(projection):
    m = re.search(r"-mRNA-(\d+)", projection)
    return int(m.group(1)) if m else 99


def read_loss(path, genes):
    status = {}
    with opener(path) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) >= 3 and t[0] == "PROJECTION" and gene_of(t[1]) in genes:
                status[t[1]] = t[2]
    return status


def read_orth(path, genes):
    """projection -> list of (ref_gene, q_gene, class); ref genes TOGA2 called one2zero."""
    orth, zero = defaultdict(list), defaultdict(set)
    with opener(path) as fh:
        next(fh, None)
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) < 5 or t[0] not in genes:
                continue
            if t[3] == "None":
                zero[t[0]].add(t[4])
            else:
                orth[t[3]].append((t[0], t[2], t[4]))
    return orth, zero


def read_bed(path, genes):
    out = {}
    with opener(path) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) >= 12 and gene_of(t[3]) in genes:
                out[t[3]] = (t[0], int(t[1]), int(t[2]), t[5], int(t[9]))
    return out


def read_fasta(path, wanted, label=""):
    """{projection: seq}; same rule as 2a for files holding a reference and a query record per id."""
    recs, cur = {}, None
    if not os.path.exists(path):
        print(f"[warn] missing {path}", file=sys.stderr)
        return {}
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                head = line[1:].strip()
                key = re.split(r"[\s|]", head)[0]
                if key in wanted:
                    cur = [head, []]
                    recs.setdefault(key, []).append(cur)
                else:
                    cur = None
            elif cur is not None:
                cur[1].append(line.strip())
    out = {}
    for key, rs in recs.items():
        if len(rs) > 1:
            q = [r for r in rs if re.search(r"query", r[0], re.I) and not re.search(r"ref", r[0], re.I)]
            if len(q) != 1:
                print(f"[warn] {label} {key}: {len(rs)} records, cannot tell which is the query; skipped",
                      file=sys.stderr)
                continue
            rs = q
        out[key] = re.sub(r"[^ACGTN]", "", "".join(rs[0][1]).upper())
    return out


def read_gene_table(path):
    with open(path) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")]
    return {r["gene_id"]: (r.get("label") or r["gene_id"]) for r in rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="e.g. ${TOGADIR}/target_HelSus")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--gene-table", help="output of 4a (gene_id and label columns)")
    g.add_argument("--genes", help="comma-separated reference gene IDs (labels = IDs)")
    ap.add_argument("--status", default="FI,I,PI", help="TOGA2 statuses allowed into the tree")
    ap.add_argument("--min-cds", type=int, default=150)
    ap.add_argument("--include-paralogs", action="store_true",
                    help="also keep projections TOGA2 did not classify as orthologs")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    label = read_gene_table(a.gene_table) if a.gene_table else {x: x for x in a.genes.split(",")}
    genes, allowed = set(label), set(a.status.split(","))
    queries = sorted(d for d in os.listdir(a.run) if d.startswith("query_"))
    os.makedirs(a.outdir, exist_ok=True)

    P, zero_rows = [], []
    for qd in queries:
        q, p = qd[6:], os.path.join(a.run, qd)
        status = read_loss(os.path.join(p, "loss_summary.tsv"), genes)
        orth, zero = read_orth(os.path.join(p, "orthology_classification.tsv"), genes)
        bed = read_bed(os.path.join(p, "query_annotation.bed"), genes)
        for gname, cls in zero.items():
            zero_rows.append((q, gname, ",".join(sorted(cls))))
        prs = sorted(set(status) | set(orth) | set(bed))
        seqs = read_fasta(os.path.join(p, "nucleotide.fa.gz"), set(prs), label=q)
        for pr in prs:
            o = orth.get(pr, [])
            chrom, s, e, strand, nex = bed.get(pr, ("NA", -1, -1, ".", 0))
            P.append(dict(query=q, ref_gene=gene_of(pr), projection=pr, chain=pr.split("#", 1)[-1],
                          status=status.get(pr, "NA"),
                          orthology_class=",".join(sorted({x[2] for x in o})) or "not_classified",
                          q_gene=",".join(sorted({x[1] for x in o})) or "NA",
                          q_chrom=chrom, q_start=s, q_end=e, q_strand=strand, n_exons=nex,
                          nt=seqs.get(pr, "")))

    pcols = ["query", "ref_gene", "projection", "chain", "status", "orthology_class", "q_gene",
             "q_chrom", "q_start", "q_end", "q_strand", "n_exons"]
    with open(os.path.join(a.outdir, "projections.tsv"), "w") as fh:
        fh.write("\t".join(pcols + ["cds_bp"]) + "\n")
        for r in P:
            fh.write("\t".join(str(r[c]) for c in pcols) + f"\t{len(r['nt'])}\n")
        for q, gname, cls in zero_rows:
            fh.write("\t".join([q, gname, "NA", "NA", "NA", cls, "None"] + ["NA"] * 5 + ["0"]) + "\n")

    # ---- units: TOGA2 ortholog genes, then (optionally) paralog loci ----
    units = []
    by_qg = defaultdict(list)
    for r in P:
        if r["q_gene"] != "NA":
            for qg in r["q_gene"].split(","):
                by_qg[(r["query"], qg)].append(r)
    for (q, qg), rs in sorted(by_qg.items()):
        units.append(dict(query=q, q_gene=qg, kind="ortholog", rs=rs))

    if a.include_paralogs:
        span = defaultdict(list)            # ortholog gene spans, to drop overlapping paralog loci
        for u in units:
            for r in u["rs"]:
                if r["q_chrom"] != "NA":
                    span[(u["query"], r["q_chrom"], r["q_strand"])].append((r["q_start"], r["q_end"]))
        para = defaultdict(list)
        for r in P:
            if r["q_gene"] == "NA" and r["q_chrom"] != "NA":
                para[(r["query"], r["q_chrom"], r["q_strand"])].append(r)
        for key, rs in sorted(para.items()):
            rs.sort(key=lambda r: r["q_start"])
            clusters, cur, end = [], [], -1
            for r in rs:
                if cur and r["q_start"] > end:
                    clusters.append(cur); cur, end = [], -1
                cur.append(r); end = max(end, r["q_end"])
            if cur:
                clusters.append(cur)
            for c in clusters:
                s, e = min(r["q_start"] for r in c), max(r["q_end"] for r in c)
                if any(os_ < e and oe > s for os_, oe in span.get(key, [])):
                    continue
                units.append(dict(query=key[0], q_gene=f"paralog:{key[1]}:{s}-{e}", kind="paralog", rs=c))

    # ---- one projection per unit, name, QC ----
    tips, used = [], set()
    for u in sorted(units, key=lambda u: (u["query"], u["kind"] != "ortholog",
                                          min(r["q_chrom"] for r in u["rs"]),
                                          min(r["q_start"] for r in u["rs"]))):
        rs = u["rs"]
        refs = sorted({r["ref_gene"] for r in rs})
        best = sorted(rs, key=lambda r: (not r["nt"], RANK.get(r["status"], 99),
                                         isoform(r["projection"]), -len(r["nt"])))[0]
        nt = best["nt"]
        aa = translate(nt)
        if aa.endswith("*"):
            nt, aa = nt[:-3], aa[:-1]
        labs = "_".join(dict.fromkeys(label[x] for x in refs))
        base = f"{u['query']}_{labs}" + ("p" if u["kind"] == "paralog" else "")
        name, k = base, 0
        while name in used:
            name = base + "bcdefghijklmnopqrstuvwxyz"[k]; k += 1
        used.add(name)
        why = []
        if not nt:
            why.append("no_sequence")
        if best["status"] not in allowed:
            why.append(f"status_{best['status']}")
        if len(nt) % 3:
            why.append("cds_not_mult3")
        if len(nt) < a.min_cds:
            why.append(f"cds<{a.min_cds}")
        tips.append(dict(tip=name, query=u["query"], kind=u["kind"], q_gene=u["q_gene"],
                         ref_genes=",".join(refs), group=labs,
                         orthology_class=",".join(sorted({r["orthology_class"] for r in rs})),
                         projection=best["projection"], status=best["status"],
                         q_chrom=best["q_chrom"], q_start=best["q_start"], q_end=best["q_end"],
                         q_strand=best["q_strand"], n_exons=best["n_exons"], cds_bp=len(nt),
                         aa_len=len(aa), internal_stops=aa.count("*"),
                         included="yes" if not why else "no", reason=";".join(why),
                         protein=aa, _nt=nt))

    cols = ["tip", "query", "kind", "q_gene", "ref_genes", "group", "orthology_class", "projection",
            "status", "q_chrom", "q_start", "q_end", "q_strand", "n_exons", "cds_bp", "aa_len",
            "internal_stops", "included", "reason", "protein"]
    with open(os.path.join(a.outdir, "tips.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for t in tips:
            fh.write("\t".join(str(t[c]) for c in cols) + "\n")
    keep = [t for t in tips if t["included"] == "yes"]
    with open(os.path.join(a.outdir, "exendin_cds.fa"), "w") as fn, \
         open(os.path.join(a.outdir, "exendin_aa.fa"), "w") as fa, \
         open(os.path.join(a.outdir, "groups.tsv"), "w") as fg:
        fg.write("tip\tgroup\n")
        for t in keep:
            fn.write(f">{t['tip']}\n{t['_nt']}\n")
            fa.write(f">{t['tip']}\n{t['protein'].replace('*', 'X')}\n")
            fg.write(f"{t['tip']}\t{t['group']}\n")

    nq = len({t["query"] for t in keep})
    print(f"[info] {len(P)} projections of {len(genes)} reference genes in {len(queries)} queries; "
          f"{len(tips)} genes/loci, {len(keep)} in the tree from {nq} species", file=sys.stderr)
    for gname in sorted(genes):
        nq_g = len({t["query"] for t in keep if gname in t["ref_genes"].split(",")})
        print(f"[info]   {gname} ({label[gname]}): orthologs in the tree from {nq_g} species", file=sys.stderr)


if __name__ == "__main__":
    main()
