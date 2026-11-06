#!/usr/bin/env python3
"""
Turn each genome's forward tblastn hits (4b) into loci, and locate the boundary genes.

Exendin-family hits are split by contig and handed to lib_scan_interval.py (the old
repo's 2f_scan_interval.py, unchanged) with --blast-out, so loci are called exactly as in
the old copy-number survey: hits clustered by position, two hits from one query kept in
one gene only if they are collinear (so tandem copies stay separate), each locus
classified as candidate_gene / N_term_only / domain_only, with premature stops and the
first residue of the hormone domain checked.

On top of 2f's tables this adds, per locus, a forward score for every query: the sum of
bitscores of that query's non-overlapping HSPs inside the locus. 4d uses it to find,
for each name group, the best forward locus in each genome (the forward half of a
reciprocal best hit).

Boundary genes (SLC4A8, AQP5) are long multi-exon genes and are not given to 2f. Their
HSPs are clustered per contig and strand (gap <= --boundary-gap-kb) and scored the same
way; every cluster scoring at least --boundary-min-frac of the best one is kept, and 4d
picks the one nearest the exendin loci (paralogues such as SLC4A10 also hit).

Usage
  python3 4c_call_loci.py --species-table species.tsv --forward-dir forward \
      --queries queries/forward_queries.faa --query-table queries/query_table.tsv \
      --boundary-queries queries/boundary_queries.faa --outdir loci

Writes to --outdir:
  forward_loci.tsv     one row per locus (all species): 2f's columns plus locus_id,
                       species, contig, best_bits, min_evalue
  forward_scores.tsv   species locus_id query group bits min_evalue
  boundary_hits.tsv    species boundary contig start end strand bits frac_of_top query_cov
  scan/<code>/         2f's per-contig input and output, kept for checking
"""
import argparse, csv, os, re, subprocess, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

SCAN = os.path.join(HERE, "lib_scan_interval.py")


def greedy_sum(hsps):
    """Sum bitscores of HSPs that do not overlap in query coordinates (best first)."""
    taken, tot, ev = [], 0.0, None
    for h in sorted(hsps, key=lambda h: -h["bits"]):
        if any(h["qs"] <= e and h["qe"] >= s for s, e in taken):
            continue
        taken.append((h["qs"], h["qe"])); tot += h["bits"]
        ev = h["ev"] if ev is None else min(ev, h["ev"])
    cov = sum(e - s + 1 for s, e in taken)
    return tot, ev, cov


def clean_contig(s):
    # older BLAST versions may report 'ref|NC_1.1|'; keep the accession
    m = re.match(r"^(?:[a-z]{2,3}\|)([^|]+)\|?$", s)
    return m.group(1) if m else s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--forward-dir", required=True, help="4b output (<code>.tblastn.tsv)")
    ap.add_argument("--queries", required=True, help="forward_queries.faa from 4a")
    ap.add_argument("--query-table", required=True, help="query_table.tsv from 4a")
    ap.add_argument("--boundary-queries", required=True, help="boundary_queries.faa from 4a")
    ap.add_argument("--min-identity", type=float, default=35.0, help="passed to 2f")
    ap.add_argument("--evalue-max", type=float, default=1e-3, help="HSPs above this are ignored")
    ap.add_argument("--cluster-kb", type=float, default=10.0, help="passed to 2f")
    ap.add_argument("--boundary-gap-kb", type=float, default=150.0)
    ap.add_argument("--boundary-min-frac", type=float, default=0.3)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    qgroup = {}
    with open(a.query_table) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            qgroup[r["id"]] = r["group"]
    blen = {n: len(s) for n, (h, s) in N.read_fasta(a.boundary_queries).items()}
    os.makedirs(a.outdir, exist_ok=True)

    L_rows, S_rows, B_rows = [], [], []
    for sp in N.read_species_table(a.species_table):
        code = sp["code"]
        path = os.path.join(a.forward_dir, f"{code}.tblastn.tsv")
        if not os.path.exists(path):
            N.info(f"{code}: no {path}; skipped")
            continue
        fam, bnd = defaultdict(list), defaultdict(list)
        n_raw = 0
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if len(f) < 13:
                continue
            n_raw += 1
            ctg = clean_contig(f[0])
            h = dict(q=f[1], qs=int(f[2]), qe=int(f[3]), ss=int(f[4]), se=int(f[5]), pid=float(f[6]),
                     ev=float(f[8]), bits=float(f[9]), raw="\t".join(f[1:13]))
            h["lo"], h["hi"], h["strand"] = min(h["ss"], h["se"]), max(h["ss"], h["se"]), \
                ("+" if h["se"] > h["ss"] else "-")
            if h["ev"] > a.evalue_max:
                continue
            if h["q"].startswith("boundary_"):
                bnd[(h["q"], ctg)].append(h)
            elif h["q"] in qgroup:
                fam[ctg].append(h)
            else:
                N.warn(f"{code}: query {h['q']} is in neither the query table nor the boundary set")

        # ---- exendin loci, via 2f per contig ----
        sdir = os.path.join(a.outdir, "scan", code)
        os.makedirs(sdir, exist_ok=True)
        nloc = 0
        for ctg in sorted(fam):
            safe = re.sub(r"[^A-Za-z0-9._-]", "_", ctg)
            bo, pre = os.path.join(sdir, safe + ".blast.tsv"), os.path.join(sdir, safe)
            with open(bo, "w") as fh:
                for h in fam[ctg]:
                    fh.write(h["raw"] + "\n")
            r = subprocess.run([sys.executable, SCAN, "--query", a.queries, "--blast-out", bo,
                                "--out-prefix", pre, "--min-identity", str(a.min_identity),
                                "--cluster-kb", str(a.cluster_kb)], capture_output=True, text=True)
            open(pre + ".log", "w").write(r.stderr)
            if r.returncode:
                N.warn(f"{code} {ctg}: lib_scan_interval.py failed:\n{r.stderr.strip()[-400:]}")
                continue
            with open(pre + "_loci.tsv") as fh:
                for L in csv.DictReader(fh, delimiter="\t"):
                    lo, hi, st = int(L["g_start"]), int(L["g_end"]), L["strand"]
                    inside = [h for h in fam[ctg] if h["strand"] == st and h["hi"] >= lo and h["lo"] <= hi]
                    byq = defaultdict(list)
                    for h in inside:
                        byq[h["q"]].append(h)
                    scores = {q: greedy_sum(hs) for q, hs in byq.items()}
                    nloc += 1
                    lid = f"{code}_L{nloc:03d}"
                    best = max(scores.items(), key=lambda kv: kv[1][0]) if scores else (None, (0, None, 0))
                    L_rows.append(dict(L, locus_id=lid, species=code, contig=ctg,
                                       best_bits=f"{best[1][0]:.1f}",
                                       min_evalue=f"{min(s[1] for s in scores.values()):.2g}" if scores else "NA"))
                    for q, (b, ev, cov) in scores.items():
                        S_rows.append([code, lid, ctg, q, qgroup[q], f"{b:.1f}", f"{ev:.2g}"])

        # ---- boundary genes ----
        nb = 0
        per_name = defaultdict(list)
        for (q, ctg), hs in bnd.items():
            name = q[len("boundary_"):].split("__")[0]
            for st in "+-":
                sh = sorted((h for h in hs if h["strand"] == st), key=lambda h: h["lo"])
                cl, end = [], -1
                for h in sh:
                    if cl and h["lo"] - end > a.boundary_gap_kb * 1000:
                        per_name[name].append((q, ctg, st, cl)); cl, end = [], -1
                    cl.append(h); end = max(end, h["hi"])
                if cl:
                    per_name[name].append((q, ctg, st, cl))
        for name, cls in per_name.items():
            scored = []
            for q, ctg, st, cl in cls:
                b, ev, cov = greedy_sum(cl)
                scored.append((b, q, ctg, st, min(h["lo"] for h in cl), max(h["hi"] for h in cl), cov))
            top = max(s[0] for s in scored)
            for b, q, ctg, st, lo, hi, cov in sorted(scored, reverse=True):
                if b >= a.boundary_min_frac * top:
                    nb += 1
                    B_rows.append([code, name, ctg, lo, hi, st, f"{b:.1f}", f"{b / top:.2f}",
                                   f"{cov / blen.get(q, cov or 1):.2f}", q])
        N.info(f"{code}: {n_raw} HSPs; {nloc} exendin-family loci on {len(fam)} contigs; "
               f"{nb} boundary cluster(s) kept for {','.join(sorted(per_name)) or 'none'}")

    lcols = ["locus_id", "species", "contig", "g_start", "g_end", "strand", "n_hits", "parts", "state",
             "best_query", "best_pident", "stops", "stops_any_hit", "domain_start", "best_bits", "min_evalue"]
    with open(os.path.join(a.outdir, "forward_loci.tsv"), "w") as fh:
        fh.write("\t".join(lcols) + "\n")
        for r in L_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in lcols) + "\n")
    with open(os.path.join(a.outdir, "forward_scores.tsv"), "w") as fh:
        fh.write("species\tlocus_id\tcontig\tquery\tgroup\tbits\tmin_evalue\n")
        for r in S_rows:
            fh.write("\t".join(r) + "\n")
    with open(os.path.join(a.outdir, "boundary_hits.tsv"), "w") as fh:
        fh.write("species\tboundary\tcontig\tstart\tend\tstrand\tbits\tfrac_of_top\tquery_cov\tquery\n")
        for r in B_rows:
            fh.write("\t".join(map(str, r)) + "\n")
    N.info(f"{len(L_rows)} loci, {len(B_rows)} boundary clusters written to {a.outdir}")


if __name__ == "__main__":
    main()
