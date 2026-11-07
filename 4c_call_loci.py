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

Split genes. 2f can call one gene as two loci, an upstream one with the signal-peptide exon
and a downstream one starting at the hormone domain (HelSus/HelCha GCGL: two pieces 712 bp
apart, the first exon modelled with intron read-through and a false stop). Neighbouring
loci on one strand within --merge-kb are joined when the queries hitting both cover
consecutive parts of themselves (the downstream part starts and ends >= --merge-min-offset
residues later in the query); tandem copies, which cover the same part twice, stay separate. Each join and each refusal is logged. The joined locus keeps the downstream
part's 2f calls (it holds the hormone domain) and gets ';merged' in its state.

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
from lib_exendin_names import info  # noqa: E402

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


def merge_partial(loci, max_gap, min_offset, code):
    """Join loci that are one gene split by 2f at an intron.

    Two loci are joined when they are neighbours on one contig and strand, at most max_gap
    apart, and the queries that hit both cover CONSECUTIVE parts of themselves: in gene
    orientation the downstream locus both starts and ends at least `min_offset` query
    residues later than the upstream one. Overlap is allowed (the upstream piece's hit can
    run on into the intron); what matters is the order. Tandem copies cover the SAME stretch
    of a query twice (same start), and a fragment beside a complete copy shares its start or
    its end, so neither is joined. Decided by majority over the queries that hit both loci.
    Returns (loci, number of joins).
    """
    def qrange(L, q):
        hs = [h for h in L["_hits"] if h["q"] == q]
        return min(h["qs"] for h in hs), max(h["qe"] for h in hs)
    joins, changed = 0, True
    while changed:
        changed = False
        loci.sort(key=lambda L: (L["_contig"], L["strand"], int(L["g_start"])))
        for i in range(len(loci) - 1):
            x, y = loci[i], loci[i + 1]
            if x["_contig"] != y["_contig"] or x["strand"] != y["strand"]:
                continue
            gap = int(y["g_start"]) - int(x["g_end"])
            if gap > max_gap:
                continue
            up, down = (x, y) if x["strand"] == "+" else (y, x)
            shared = {h["q"] for h in up["_hits"]} & {h["q"] for h in down["_hits"]}
            if not shared:
                continue
            consec = []
            for q in shared:
                (us, ue), (ds, de) = qrange(up, q), qrange(down, q)
                consec.append(ds - us >= min_offset and de - ue >= min_offset)
            ok = sum(consec) > len(consec) / 2
            ids = f"{x['_contig']}:{x['g_start']}-{x['g_end']} + {y['g_start']}-{y['g_end']} ({x['strand']}, gap {gap} bp)"
            if not ok:
                info(f"{code}: not joined, they cover the same part of {len(consec) - sum(consec)}/{len(consec)} "
                     f"shared queries (tandem copies): {ids}")
                continue
            m = dict(down)                       # the part with the hormone domain keeps its 2f calls
            m["g_start"] = str(min(int(x["g_start"]), int(y["g_start"])))
            m["g_end"] = str(max(int(x["g_end"]), int(y["g_end"])))
            m["n_hits"] = str(int(x.get("n_hits", 0) or 0) + int(y.get("n_hits", 0) or 0))
            m["parts"] = f"{up.get('parts', '')}|{down.get('parts', '')}"
            st = down["state"].split(";")
            m["state"] = ";".join(["candidate_gene"] + [f for f in st[1:] if f != "merged"] + ["merged"])
            m["_hits"] = up["_hits"] + down["_hits"]
            loci[i:i + 2] = [m]
            joins += 1
            info(f"{code}: joined, consecutive parts of {sum(consec)}/{len(consec)} shared queries "
                 f"(one gene split at an intron): {ids}")
            changed = True
            break
    return loci, joins


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
    ap.add_argument("--merge-kb", type=float, default=5.0,
                    help="join neighbouring loci on one strand at most this far apart when they cover "
                         "consecutive parts of the same queries (one gene split at an intron); 0 = never")
    ap.add_argument("--merge-min-offset", type=int, default=10,
                    help="the downstream part must start and end this many query residues later than the "
                         "upstream part (consecutive parts of one gene, not two copies of the same part)")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    qgroup = {}
    with open(a.query_table) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            qgroup[r["id"]] = r["group"]
    blen = {n: len(s) for n, (h, s) in N.read_fasta(a.boundary_queries).items()}
    os.makedirs(a.outdir, exist_ok=True)

    L_rows, S_rows, B_rows = [], [], []
    nmerged = 0
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
        nloc, sp_loci = 0, []
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
                    L["_hits"] = [h for h in fam[ctg] if h["strand"] == st and h["hi"] >= lo and h["lo"] <= hi]
                    L["_contig"] = ctg
                    sp_loci.append(L)

        # ---- join partial loci of one gene (2f can split a gene at an intron) ----
        if a.merge_kb > 0:
            sp_loci, merged = merge_partial(sp_loci, a.merge_kb * 1000, a.merge_min_offset, code)
            nmerged += merged
        sp_loci.sort(key=lambda L: (L["_contig"], int(L["g_start"])))
        for L in sp_loci:
            byq = defaultdict(list)
            for h in L["_hits"]:
                byq[h["q"]].append(h)
            scores = {q: greedy_sum(hs) for q, hs in byq.items()}
            nloc += 1
            lid = f"{code}_L{nloc:03d}"
            best = max(scores.items(), key=lambda kv: kv[1][0]) if scores else (None, (0, None, 0))
            row = {k: v for k, v in L.items() if not k.startswith("_")}
            if best[0]:
                row["best_query"] = best[0]
            L_rows.append(dict(row, locus_id=lid, species=code, contig=L["_contig"],
                               best_bits=f"{best[1][0]:.1f}",
                               min_evalue=f"{min(s[1] for s in scores.values()):.2g}" if scores else "NA"))
            for q, (b, ev, cov) in scores.items():
                S_rows.append([code, lid, L["_contig"], q, qgroup[q], f"{b:.1f}", f"{ev:.2g}"])

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
    N.info(f"{len(L_rows)} loci ({nmerged} joined from split genes), {len(B_rows)} boundary clusters "
           f"written to {a.outdir}")


if __name__ == "__main__":
    main()
