#!/usr/bin/env python3
"""
Protein and CDS for every curated locus (4d curated_loci.tsv), for gene trees later.

Each locus is rebuilt by lib_rescue_locus.py (the old repo's 2e_rescue_locus.py,
unchanged) in protein mode: the locus's best forward query of its own name group is
aligned to the locus (+/- --flank) with tblastn, collinear hits are spliced, and the
longest ORF is reported. The hormone domain is then located with 2f's motif, and its
first residue reported (every functional peptide in this family starts His or Tyr).

Models are evidence-based but automatic: check the exon table of anything whose length
is far from 84-113 aa, as in the old step-2 QC.

Usage
  python3 4e_build_models.py --species-table species.tsv --genomedir ${GENOMEDIR} \
      --curated reciprocal/curated_loci.tsv --queries queries/forward_queries.faa \
      --outdir models

Writes to --outdir: <name>.{aa,cds}.fa and <name>.exons.tsv from 2e, plus
  curated_aa.fa, curated_cds.fa   all models, named as in curated_loci.tsv
  models.tsv                      name, query used, exons, aa length, domain start, status
"""
import argparse, csv, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

RESCUE = os.path.join(HERE, "lib_rescue_locus.py")
MOTIF = re.compile(r"[HY][ASG][DEQ][AGP][ATVIGSM]F[TSI]")      # as in 2f / 2d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--genomedir", required=True)
    ap.add_argument("--curated", required=True, help="curated_loci.tsv from 4d")
    ap.add_argument("--queries", required=True, help="forward_queries.faa from 4a")
    ap.add_argument("--flank", type=int, default=2000)
    ap.add_argument("--min-identity", type=float, default=30.0)
    ap.add_argument("--evalue", default="1e-2")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    table = {r["code"]: r for r in N.read_species_table(a.species_table)}
    qs = N.read_fasta(a.queries)
    with open(a.curated) as fh:
        cur = list(csv.DictReader(fh, delimiter="\t"))
    rows, aa_out, cds_out = [], [], []
    fai_cache = {}
    for r in cur:
        code, name, q = r["species"], r["name"], r["best_forward_query"]
        files = N.find_genome_files(a.genomedir, code, table[code]["accession"])
        if files["fasta"] not in fai_cache:
            fai_cache[files["fasta"]] = {l.split("\t")[0]: int(l.split("\t")[1])
                                         for l in open(files["fasta"] + ".fai")}
        clen = fai_cache[files["fasta"]].get(r["contig"], 10 ** 12)
        lo, hi = max(1, int(r["g_start"]) - a.flank), min(clen, int(r["g_end"]) + a.flank)
        if q not in qs:
            rows.append([name, q, 0, 0, "", "query_missing"]); continue
        qf = os.path.join(a.outdir, f"{name}.query.faa")
        open(qf, "w").write(f">{q}\n{qs[q][1]}\n")
        p = subprocess.run([sys.executable, RESCUE, "--genome", files["fasta"], "--contig", r["contig"],
                            "--region", f"{lo}-{hi}", "--query", qf, "--protein", "--name", name,
                            "--outdir", a.outdir, "--min-identity", str(a.min_identity),
                            "--evalue", a.evalue], capture_output=True, text=True)
        open(os.path.join(a.outdir, f"{name}.log"), "w").write(p.stderr)
        if p.returncode:
            msg = (p.stderr.strip().splitlines() or ["failed"])[-1]
            rows.append([name, q, 0, 0, "", "failed: " + msg[:80]])
            N.warn(f"{name}: model failed: {msg}")
            continue
        aa = N.read_fasta(os.path.join(a.outdir, f"{name}.aa.fa"))[name][1]
        cds = N.read_fasta(os.path.join(a.outdir, f"{name}.cds.fa"))[name][1]
        nex = sum(1 for _ in open(os.path.join(a.outdir, f"{name}.exons.tsv"))) - 1
        m = MOTIF.search(aa)
        dom = f"{aa[m.start()]}@{m.start() + 1}" if m else "none"
        status = "ok" if m else "no_hormone_motif"
        if not (60 <= len(aa) <= 150):
            status += f";length_{len(aa)}aa"
        rows.append([name, q, nex, len(aa), dom, status])
        aa_out.append((name, aa)); cds_out.append((name, cds))
    with open(os.path.join(a.outdir, "models.tsv"), "w") as fh:
        fh.write("name\tquery\texons\taa_len\tdomain_start\tstatus\n")
        for row in rows:
            fh.write("\t".join(map(str, row)) + "\n")
    for fn, recs in (("curated_aa.fa", aa_out), ("curated_cds.fa", cds_out)):
        with open(os.path.join(a.outdir, fn), "w") as fh:
            for n, s in recs:
                fh.write(f">{n}\n{s}\n")
    N.info(f"{len(aa_out)}/{len(cur)} loci modelled; "
           f"{sum(1 for r in rows if r[5] != 'ok')} need a look (models.tsv status != ok)")


if __name__ == "__main__":
    main()
