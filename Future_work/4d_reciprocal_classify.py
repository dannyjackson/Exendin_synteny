#!/usr/bin/env python3
"""
Reverse half of the reciprocal BLAST, and the classification of every forward locus
(4c) into orthologs, paralogs and rejects.

Reverse search. Each locus (+/- --flank bp) is cut from its genome with samtools and
searched with blastx against the reverse database from 4a (whole reference proteome(s)
plus all exendin-family queries). For each locus and database protein the score is the
sum of bitscores of non-overlapping HSPs on the locus's strand, so a three-exon gene is
scored as a whole. Query-set proteins from the locus's own species are ignored, so a
locus cannot be named simply after its own annotation; reference-proteome proteins are
always used, as they are the competitors (proglucagon, GIP, VIP, PACAP ...).

Classification of a locus in species X
  rejected        best reverse hit is not an exendin-family protein (e.g. proglucagon)
  no_reverse_hit  nothing in the reverse database hits it
  ortholog        best reverse hit is family group G, AND this locus is X's best forward
                  locus for the queries of group G: a reciprocal best hit for G
  paralog         best reverse hit is family group G, but another locus of X is the best
                  forward hit for G (an extra copy named after G)
Confidence (ortholog and paralog only)
  high    margin (best family reverse score / best non-family reverse score) >= --min-margin,
          forward e-value <= --max-evalue, and 2f state 'candidate_gene' (both the exon-1
          and the hormone-domain part of the precursor found)
  medium  any of those fails (remnants, weak hits, ambiguous reverse hit)
'disabled' is reported separately (2f premature_stop / domain_starts_X): a disabled copy
can still be a high-confidence paralog; it is drawn dashed in 4f.

Placement, as in 2a: interval / near_boundary / unplaced / elsewhere, using the boundary
clusters of 4c. A boundary is picked as in 2a: the SLC4A8-AQP5 pair on one contig with
the most accepted loci between them, else each one nearest accepted loci.

Usage
  python3 4d_reciprocal_classify.py --species-table species.tsv --genomedir ${GENOMEDIR} \
      --loci-dir loci --queries-dir queries --threads 8 --outdir reciprocal

  # re-classify without re-running blastx (e.g. with another --min-margin)
  python3 4d_reciprocal_classify.py ... --reuse-blastx

Writes to --outdir:
  reverse_db.*            BLAST protein database
  regions/<code>.fa       locus sequences; regions/<code>.blastx.tsv  reverse hits
  reciprocal_loci.tsv     every forward locus with its reverse hits, class, confidence,
                          placement and name
  curated_loci.tsv        the high-confidence orthologs and paralogs (the curated set)
  boundaries_used.tsv     the boundary clusters used for placement
  summary.tsv             per species and group: orthologs, paralogs (high / medium)
"""
import argparse, csv, os, subprocess, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

COLS = ["name", "locus_id", "species", "contig", "g_start", "g_end", "strand", "group", "class", "confidence",
        "disabled", "placement", "state", "parts", "domain_start", "best_forward_query", "forward_bits",
        "forward_min_evalue", "forward_best_locus_for_group", "rev_best", "rev_best_group", "rev_best_bits",
        "rev_best_desc", "rev_best_nonfamily", "rev_best_nonfamily_desc", "rev_nonfamily_bits", "margin",
        "reason"]


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode:
        sys.exit(f"[error] {cmd[0]} failed:\n{r.stderr.strip()[-600:]}")
    return r.stdout


def greedy_sum(hsps):
    taken, tot = [], 0.0
    for s, e, b in sorted(hsps, key=lambda h: -h[2]):
        if any(s <= te and e >= ts for ts, te in taken):
            continue
        taken.append((s, e)); tot += b
    return tot


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--genomedir", required=True)
    ap.add_argument("--loci-dir", required=True, help="4c output")
    ap.add_argument("--queries-dir", required=True, help="4a output")
    ap.add_argument("--flank", type=int, default=1000, help="bp added to each side of a locus for blastx")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--evalue", default="1e-3", help="blastx e-value cut-off")
    ap.add_argument("--min-margin", type=float, default=1.2)
    ap.add_argument("--max-evalue", type=float, default=1e-5)
    ap.add_argument("--boundary-frac", type=float, default=0.5,
                    help="boundary clusters scoring below this fraction of the best are not used")
    ap.add_argument("--pad-kb", type=float, default=500.0)
    ap.add_argument("--max-dist-kb", type=float, default=1000.0)
    ap.add_argument("--max-interval-kb", type=float, default=5000.0)
    ap.add_argument("--reuse-blastx", action="store_true", help="use regions/<code>.blastx.tsv if present")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    os.makedirs(os.path.join(a.outdir, "regions"), exist_ok=True)
    rev = {r["id"]: r for r in read_tsv(os.path.join(a.queries_dir, "reverse_table.tsv"))}
    qgroup = {r["id"]: r["group"] for r in read_tsv(os.path.join(a.queries_dir, "query_table.tsv"))}
    loci = read_tsv(os.path.join(a.loci_dir, "forward_loci.tsv"))
    scores = read_tsv(os.path.join(a.loci_dir, "forward_scores.tsv"))
    bhits = read_tsv(os.path.join(a.loci_dir, "boundary_hits.tsv"))

    db = os.path.join(a.outdir, "reverse_db")
    if not (a.reuse_blastx and all(os.path.exists(os.path.join(a.outdir, "regions", f"{c}.blastx.tsv"))
                                   for c in {L["species"] for L in loci})):
        if not os.path.exists(db + ".pdb") and not os.path.exists(db + ".pin"):
            run(["makeblastdb", "-in", os.path.join(a.queries_dir, "reverse_db.faa"), "-dbtype", "prot",
                 "-parse_seqids", "-out", db])

    by_sp = defaultdict(list)
    for L in loci:
        by_sp[L["species"]].append(L)
    table = {r["code"]: r for r in N.read_species_table(a.species_table)}

    # ---------------- reverse search ----------------
    rev_scores = {}            # locus_id -> {subject: score}
    for code, Ls in sorted(by_sp.items()):
        bx = os.path.join(a.outdir, "regions", f"{code}.blastx.tsv")
        if not (a.reuse_blastx and os.path.exists(bx)):
            files = N.find_genome_files(a.genomedir, code, table[code]["accession"])
            fai = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(files["fasta"] + ".fai")}
            regs = []
            for L in Ls:
                lo = max(1, int(L["g_start"]) - a.flank)
                hi = min(fai.get(L["contig"], int(L["g_end"]) + a.flank), int(L["g_end"]) + a.flank)
                regs.append(f"{L['contig']}:{lo}-{hi}")
                L["_reg_lo"] = lo
            rfile = os.path.join(a.outdir, "regions", f"{code}.regions.txt")
            open(rfile, "w").write("\n".join(regs) + "\n")
            raw = run(["samtools", "faidx", files["fasta"], "-r", rfile])
            fa, k = [], -1
            for line in raw.splitlines():
                if line.startswith(">"):
                    k += 1; fa.append(f">{Ls[k]['locus_id']}")
                else:
                    fa.append(line)
            qfa = os.path.join(a.outdir, "regions", f"{code}.fa")
            open(qfa, "w").write("\n".join(fa) + "\n")
            out = run(["blastx", "-query", qfa, "-db", db, "-evalue", a.evalue, "-seg", "no",
                       "-max_target_seqs", "200", "-num_threads", str(a.threads), "-outfmt",
                       "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"])
            open(bx, "w").write(out)
        strand = {L["locus_id"]: L["strand"] for L in Ls}
        hs = defaultdict(lambda: defaultdict(list))
        for line in open(bx):
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            lid, sub, qs, qe, bits = f[0], f[1], int(f[4]), int(f[5]), float(f[9])
            if lid not in strand or sub not in rev:
                continue
            if ("+" if qe > qs else "-") != strand[lid]:
                continue                                 # other strand: not this gene
            r = rev[sub]
            if r["source"] != "proteome" and r["species"] == code:
                continue                                 # own species' query: not evidence
            hs[lid][sub].append((min(qs, qe), max(qs, qe), bits))
        for lid in strand:
            rev_scores[lid] = {s: greedy_sum(h) for s, h in hs[lid].items()}
        N.info(f"{code}: reverse hits for {sum(1 for l in strand if rev_scores[l])}/{len(strand)} loci")

    # ---------------- forward best locus per group ----------------
    fwd = defaultdict(dict)            # (species, group) -> {locus: best bits over the group's queries}
    best_q = {}                        # (locus, group) -> (query, bits)
    for s in scores:
        g = s["group"] if s["query"] not in qgroup else qgroup[s["query"]]
        b = float(s["bits"])
        k = (s["species"], g)
        if b > fwd[k].get(s["locus_id"], -1):
            fwd[k][s["locus_id"]] = b
        if b > best_q.get((s["locus_id"], g), ("", -1))[1]:
            best_q[(s["locus_id"], g)] = (s["query"], b)
    fwd_best = {k: max(v, key=v.get) for k, v in fwd.items()}

    # ---------------- classify ----------------
    out = []
    for L in loci:
        lid, code = L["locus_id"], L["species"]
        sc = rev_scores.get(lid, {})
        r = dict(locus_id=lid, species=code, contig=L["contig"], g_start=L["g_start"], g_end=L["g_end"],
                 strand=L["strand"], state=L["state"], parts=L["parts"], domain_start=L["domain_start"],
                 forward_bits=L["best_bits"], forward_min_evalue=L["min_evalue"],
                 disabled="yes" if ("premature_stop" in L["state"] or "domain_starts" in L["state"]) else "no")
        fam = {s: v for s, v in sc.items() if N.is_family(rev[s]["group"])}
        non = {s: v for s, v in sc.items() if not N.is_family(rev[s]["group"])}
        bf = max(fam, key=fam.get) if fam else None
        bn = max(non, key=non.get) if non else None
        top = max(sc, key=sc.get) if sc else None
        r.update(rev_best=top or "NA", rev_best_group=rev[top]["group"] if top else "NA",
                 rev_best_bits=f"{sc[top]:.1f}" if top else "0", rev_best_desc=rev[top]["description"] if top else "",
                 rev_best_nonfamily=bn or "NA", rev_best_nonfamily_desc=rev[bn]["description"] if bn else "",
                 rev_nonfamily_bits=f"{non[bn]:.1f}" if bn else "0")
        margin = (fam[bf] / non[bn]) if (bf and bn) else (float("inf") if bf else 0.0)
        r["margin"] = "inf" if margin == float("inf") else f"{margin:.2f}"
        why = []
        if not top:
            r.update(cls="no_reverse_hit", group="NA"); why.append("no reverse hit")
        elif not N.is_family(rev[top]["group"]):
            r.update(cls="rejected", group=rev[top]["group"])
            why.append(f"best reverse hit {rev[top]['group']} ({rev[top]['description'][:40]})")
        else:
            g = rev[top]["group"]
            r["group"] = g
            fb = fwd_best.get((code, g))
            r["forward_best_locus_for_group"] = fb or "none"
            r["cls"] = "ortholog" if fb == lid else "paralog"
            if fb != lid:
                why.append(f"best forward locus for {g} is {fb}" if fb else f"no forward hit from {g} queries")
        q = best_q.get((lid, r["group"]))
        r["best_forward_query"] = q[0] if q else L["best_query"]
        if r["cls"] in ("ortholog", "paralog"):
            ok = []
            if margin < a.min_margin:
                ok.append(f"margin {margin:.2f} < {a.min_margin}")
            ev = float(L["min_evalue"]) if L["min_evalue"] not in ("", "NA") else 1.0
            if ev > a.max_evalue:
                ok.append(f"forward e-value {ev:.2g} > {a.max_evalue:g}")
            if not L["state"].startswith("candidate_gene"):
                ok.append(f"2f state {L['state'].split(';')[0]}")
            r["confidence"] = "medium" if ok else "high"
            why += ok
        else:
            r["confidence"] = "none"
        r["class"] = r.pop("cls")
        r["reason"] = "; ".join(why) or "-"
        out.append(r)

    # ---------------- boundaries and placement ----------------
    pad, maxd, maxi = a.pad_kb * 1000, a.max_dist_kb * 1000, a.max_interval_kb * 1000
    used_b = []
    for code, Ls in by_sp.items():
        acc = [r for r in out if r["species"] == code and r["class"] in ("ortholog", "paralog")]
        bh = [b for b in bhits if b["species"] == code and float(b["frac_of_top"]) >= a.boundary_frac]
        names = sorted({b["boundary"] for b in bhits})
        cand = {n: [b for b in bh if b["boundary"] == n] for n in names}

        def between(c, lo, hi):
            return sum(1 for r in acc if r["contig"] == c and int(r["g_end"]) >= lo and int(r["g_start"]) <= hi)

        def near(b):
            return sum(1 for r in acc if r["contig"] == b["contig"] and
                       min(abs(int(r["g_start"]) - int(b["end"])), abs(int(b["start"]) - int(r["g_end"]))) <= maxd)
        chosen, pair = {}, False
        if len(names) == 2:
            n0, n1 = names
            pairs = []
            for s in cand[n0]:
                for t in cand[n1]:
                    if s["contig"] != t["contig"]:
                        continue
                    lo, hi = min(int(s["start"]), int(t["start"])), max(int(s["end"]), int(t["end"]))
                    if hi - lo <= maxi:
                        pairs.append(((-between(s["contig"], lo, hi), -(float(s["bits"]) + float(t["bits"])),
                                       hi - lo), s, t))
            if pairs:
                _, s, t = min(pairs, key=lambda p: p[0])
                chosen, pair = {n0: s, n1: t}, True
        if not chosen:
            for n in names:
                if cand[n]:
                    chosen[n] = min(cand[n], key=lambda b: (-near(b), -float(b["bits"])))
        for n, b in chosen.items():
            used_b.append(dict(b, placement="boundary_pair" if pair else "boundary_single"))
        bctg = defaultdict(list)
        for b in chosen.values():
            bctg[b["contig"]].append(b)
        if pair:
            pctg = next(iter(chosen.values()))["contig"]
            lo = min(int(b["start"]) for b in chosen.values()) - pad
            hi = max(int(b["end"]) for b in chosen.values()) + pad
        for r in (r for r in out if r["species"] == code):
            s, e = int(r["g_start"]), int(r["g_end"])
            if pair and r["contig"] == pctg and e >= lo and s <= hi:
                r["placement"] = "interval"
            elif r["contig"] in bctg and any(min(abs(s - int(b["end"])), abs(int(b["start"]) - e)) <= maxd
                                             for b in bctg[r["contig"]]):
                r["placement"] = "near_boundary"
            elif r["contig"] not in bctg:
                r["placement"] = "unplaced"
            else:
                r["placement"] = "elsewhere"
        N.info(f"{code}: boundaries " + (", ".join(f"{n} {b['contig']}:{b['start']}-{b['end']}"
                                                  for n, b in chosen.items()) or "not found")
               + (" (pair)" if pair else ""))

    # ---------------- names ----------------
    k = defaultdict(int)
    for r in sorted(out, key=lambda r: (r["species"], r["contig"], int(r["g_start"]))):
        if r["class"] == "ortholog":
            r["name"] = f"{r['species']}_{N.group_abbrev(r['group'])}"
        elif r["class"] == "paralog":
            k[(r["species"], r["group"])] += 1
            r["name"] = f"{r['species']}_{N.group_abbrev(r['group'])}p{k[(r['species'], r['group'])]}"
        else:
            r["name"] = r["locus_id"]

    def write(fn, rows):
        with open(os.path.join(a.outdir, fn), "w") as fh:
            fh.write("\t".join(COLS) + "\n")
            for r in rows:
                fh.write("\t".join(str(r.get(c, "")) for c in COLS) + "\n")
    out.sort(key=lambda r: (r["species"], r["contig"], int(r["g_start"])))
    write("reciprocal_loci.tsv", out)
    cur = [r for r in out if r["confidence"] == "high"]
    write("curated_loci.tsv", cur)
    with open(os.path.join(a.outdir, "boundaries_used.tsv"), "w") as fh:
        bc = ["species", "boundary", "contig", "start", "end", "strand", "bits", "frac_of_top", "query_cov",
              "query", "placement"]
        fh.write("\t".join(bc) + "\n")
        for b in used_b:
            fh.write("\t".join(str(b.get(c, "")) for c in bc) + "\n")
    summ = defaultdict(lambda: defaultdict(int))
    for r in out:
        if r["class"] in ("ortholog", "paralog"):
            summ[(r["species"], r["group"])][f"{r['class']}_{r['confidence']}"] += 1
    with open(os.path.join(a.outdir, "summary.tsv"), "w") as fh:
        ks = ["ortholog_high", "paralog_high", "ortholog_medium", "paralog_medium"]
        fh.write("species\tgroup\t" + "\t".join(ks) + "\n")
        for (sp, g), d in sorted(summ.items()):
            fh.write(f"{sp}\t{g}\t" + "\t".join(str(d[x]) for x in ks) + "\n")
    cls = defaultdict(int)
    for r in out:
        cls[f"{r['class']}/{r['confidence']}"] += 1
    N.info(f"{len(out)} loci: " + ", ".join(f"{c} {n}" for c, n in sorted(cls.items())))
    N.info(f"{len(cur)} high-confidence loci in curated_loci.tsv")


if __name__ == "__main__":
    main()
