#!/usr/bin/env python3
"""
The genome annotations behind every exendin-family gene of one genome (default HelSus):
for each locus, the genes and transcripts of each annotation that overlap it, so the
loci can be paired with expression data quantified against those annotations
(10_toga2_exendin_tree.md, step 6).

  --gff LABEL=PATH   an annotation of the same assembly, e.g.
                     EviAnn=${HELSUS_EVIANN_GFF}   (the one TOGA2 projected from)
                     NCBI=.../GCF_057233095.1/genomic.gff
                     Given as many times as needed.

Loci come from names/names.tsv (10f): coordinates of the TOGA2 projection (or of the
rebuilt model). These are in the TOGA2 query's sequence names (HelSus: 'chr2'). An
annotation that uses other names (NCBI: 'NC_139459.1') is matched through its NCBI
sequence report (sequence_report.jsonl next to the GFF: chromosome '2' <-> chr2) or
--chrom-map chr2=NC_139459.1.

For each locus and annotation: every gene (and pseudogene) overlapping it, with strand,
overlap in bp and as % of the locus, and its transcripts (transcript IDs, products).
'none' = the annotation has no gene there: the locus is unannotated in it. A gene much
longer than the locus that overlaps it on the other strand is usually a host gene, not
the exendin: the strand and overlap columns tell them apart.

Rows: GCGL, Exendin_anc and lettered copies of --species, in chromosome order. The first
column is the name from 10f, as a key for joining; every other column comes from the
annotations. Writes --out (TSV, one row per locus x annotation gene) and prints a summary.

Usage
  python3 10g_gene_table.py --species HelSus --names names/names.tsv \
      --gff EviAnn=${HELSUS_EVIANN_GFF} --gff NCBI=${HELSUS_NCBI_GFF} --out names/HelSus_annotations.tsv
"""
import argparse, csv, glob, gzip, json, os, re, sys
from collections import defaultdict
from urllib.parse import unquote


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def seq_report_map(gff):
    """{'chr2': 'NC_139459.1', '2': 'NC_139459.1', 'CM...': 'NC_...'} from an NCBI sequence report next to the GFF."""
    out = {}
    for p in glob.glob(os.path.join(os.path.dirname(gff), "sequence_report.jsonl")):
        for line in open(p):
            try:
                d = json.loads(line)
            except json.JSONDecodeError:
                continue
            acc = d.get("refseqAccession") or d.get("genbankAccession")
            for k in (d.get("chrName"), d.get("sequenceName"), d.get("genbankAccession"), d.get("ucscStyleName")):
                if k and acc:
                    out[k] = acc
                    if k == d.get("chrName"):
                        out["chr" + k] = acc
    return out


def load_genes(gff, contigs):
    """genes on the given contigs: {contig: [gene dicts with transcripts]}."""
    genes, tx = {}, defaultdict(list)
    with opener(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 9 or t[0] not in contigs:
                continue
            at = {k: unquote(v) for k, v in (kv.split("=", 1) for kv in t[8].strip().strip(";").split(";") if "=" in kv)}
            if t[2] in ("gene", "pseudogene"):
                gid = at.get("ID", "")
                dbx = dict(x.split(":", 1) for x in at.get("Dbxref", "").split(",") if ":" in x)
                genes[gid] = dict(contig=t[0], start=int(t[3]), end=int(t[4]), strand=t[6], id=gid,
                                  name=at.get("gene") or at.get("Name") or gid, geneid=dbx.get("GeneID", ""),
                                  biotype=at.get("gene_biotype", t[2]),
                                  desc=at.get("description") or at.get("Note") or at.get("product", ""))
            elif t[2] in ("mRNA", "transcript", "lnc_RNA", "ncRNA", "primary_transcript") and at.get("Parent"):
                tx[at["Parent"].split(",")[0]].append(
                    (at.get("transcript_id") or at.get("Name") or at.get("ID", ""), at.get("product") or at.get("Note", "")))
    by = defaultdict(list)
    for g in genes.values():
        g["tx"] = tx.get(g["id"], [])
        by[g["contig"]].append(g)
    return by


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--species", default="HelSus")
    ap.add_argument("--names", required=True, help="names/names.tsv from 10f")
    ap.add_argument("--gff", action="append", required=True, help="LABEL=PATH of an annotation of this assembly")
    ap.add_argument("--chrom-map", action="append", default=[], help="LOCUS_NAME=GFF_NAME, e.g. chr2=NC_139459.1")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    loci = [r for r in read_tsv(a.names) if r["species"] == a.species and r["class"] not in ("GCG", "GIP")
            and r.get("start", "").isdigit()]
    if not loci:
        sys.exit(f"[error] no {a.species} loci with coordinates in {a.names}")
    loci.sort(key=lambda r: (r["contig"], int(r["start"])))
    manual = dict(m.split("=", 1) for m in a.chrom_map)

    rows, summary = [], defaultdict(dict)
    for spec in a.gff:
        label, _, path = spec.partition("=")
        sr = seq_report_map(path)
        alias = {}
        for r in loci:
            c = r["contig"]
            c2 = re.sub(r"__(\d+)$", r".\1", c)
            alias[c] = manual.get(c) or manual.get(c2) or (c if c in sr.values() else None) or sr.get(c) or sr.get(c2) or c2
        genes = load_genes(path, set(alias.values()))
        n_on = sum(len(v) for v in genes.values())
        if not n_on:
            print(f"[warn] {label}: no genes on {', '.join(sorted(set(alias.values())))} in {path}; "
                  f"give the sequence name with --chrom-map {loci[0]['contig']}=<name in the GFF>", file=sys.stderr)
        else:
            print(f"[info] {label}: {path}; locus sequences mapped as "
                  + ", ".join(f"{k}->{v}" for k, v in sorted(alias.items()) if k != v or True), file=sys.stderr)
        for r in loci:
            s, e, ctg = int(r["start"]), int(r["end"]), alias[r["contig"]]
            hits = []
            for g in genes.get(ctg, []):
                ov = min(e, g["end"]) - max(s, g["start"]) + 1
                if ov > 0:
                    hits.append((ov, g))
            hits.sort(key=lambda h: -h[0])
            if not hits:
                rows.append([r["new_name"], r["contig"], s, e, label, "none", "", "", "", "", "", "", "", ""])
                summary[r["new_name"]][label] = "none"
                continue
            for ov, g in hits:
                rows.append([r["new_name"], r["contig"], s, e, label, g["name"], g["id"], g["geneid"], g["biotype"],
                             g["strand"], f"{g['start']}-{g['end']}", ov, f"{100 * ov / (e - s + 1):.0f}",
                             g["desc"], ";".join(t for t, _ in g["tx"]), ";".join(sorted({p for _, p in g["tx"] if p}))])
            def host(g):        # much longer than the locus: probably a gene the locus sits inside, not the exendin
                return (g["end"] - g["start"] + 1) > 5 * (e - s + 1)
            summary[r["new_name"]][label] = ", ".join(
                f"{g['name']}{' (GeneID ' + g['geneid'] + ')' if g['geneid'] else ''}"
                f"{' [' + str(len(g['tx'])) + ' tx]' if g['tx'] else ''}"
                f"{' (spans locus: host gene?)' if host(g) else ''}" for _, g in hits)

    cols = ["locus_name", "locus_contig", "locus_start", "locus_end", "annotation", "gene", "gene_id", "ncbi_geneid",
            "biotype", "gene_strand", "gene_span", "overlap_bp", "overlap_pct_of_locus", "description",
            "transcripts", "transcript_products"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            r = r + [""] * (len(cols) - len(r))
            fh.write("\t".join(map(str, r)) + "\n")
    labels = [g.partition("=")[0] for g in a.gff]
    w = max(len(r["new_name"]) for r in loci) + 2
    print("locus".ljust(w) + "position".ljust(28) + "  ".join(labels))
    for r in loci:
        pos = f"{r['contig']}:{r['start']}-{r['end']}"
        print(r["new_name"].ljust(w) + pos.ljust(28) + "  |  ".join(summary[r["new_name"]].get(l, "-") for l in labels))
    print(f"[info] {len(loci)} {a.species} loci; full table (transcript IDs, products, overlaps) in {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
