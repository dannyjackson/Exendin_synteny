#!/usr/bin/env python3
"""
Collect every exendin gene that TOGA2 annotated, in every query of every TOGA2 run, for a
curated alignment and tree (10_toga2_exendin_tree.md).

Gene set. A projection is an exendin if its REFERENCE gene is named exendin
(lib_exendin_names.name_group: exendin-N, exendin-N-like, unnumbered exendin), using each
run's reference annotation (--ref-gff RUN=GFF) or a name map. --include-gcgl adds genes
named glucagon-like (GCGL), the ancestral gene, as an outgroup for rooting.

Filters (TOGA2 output only, nothing rebuilt):
  ortholog projections   status in --status (default FI,I = intact)
  paralogous projections (#paralog, or not classified as orthologs) status in
                         --paralog-status (default: any; see the log for each one's status)
  retrogene / processed-pseudogene projections: left out
  every kept sequence: CDS present, length a multiple of 3, >= --min-cds
Internal stops are reported, not filtered; that is a curation call.

Loci. Projections of one run that overlap on one strand are one locus (isoforms, or the
same locus reached from two reference exendins); the projection with the best status, then
the lowest isoform number, then the longest CDS represents it. The same locus found by
several runs (overlap in the same query genome) is kept once: best status, then run order
of --run, then longest CDS; the other runs are listed in 'also_in'.

Tip names: <Query>_<label>[p][b..] -- label from the reference exendin(s) (Ex2, Ex4,
Ex2_4, ExN, GCGL), 'p' for a paralogous projection, letters for further copies.

Read-only on the TOGA2 directory. Writes to --outdir:
  candidates.tsv   every exendin locus found, kept or not, with the reason
  exendin_cds.fa   CDS of kept loci (final stop removed)
  exendin_aa.fa    their translations (internal stops as X)
  curation.tsv     one row per kept tip, action 'keep', with automatic notes; edit it and
                   apply it with 10b_apply_curation.py

Usage
  python3 10a_collect_toga2_exendins.py \
      --run HelSus="${TOGADIR}/target_HelSus/query_{q}" \
      --run HemCap="${TOGADIR}/target_HemCap/query_{q}" \
      --run AnoCar="${TOGADIR}/query_{q}/target_AnoCar/output" \
      --run CroPor="${TOGADIR}/query_{q}/target_CroPor/output" \
      --ref-gff HelSus=${HELSUS_EVIANN_GFF} --ref-gff HemCap=${HEMCAP_GFF} ... \
      --outdir collect 2> collect.log
"""
import argparse, glob, importlib.util, os, re, sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402


def load(fn, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


T3C = load("3c_plot_toga2_synteny.py", "toga2_interval")    # names, projections, status
T3B = load("3b_collect_toga2_exendins.py", "toga2_collect")   # FASTA with ref/query records, translate
RANK = T3B.RANK


def queries_of(pattern):
    rx = re.compile("^" + re.escape(pattern).replace(re.escape("{q}"), "([^/]+)") + "/?$")
    return sorted(m.group(1) for m in (rx.match(p.rstrip("/")) for p in glob.glob(pattern.replace("{q}", "*")))
                  if m and os.path.isdir(m.group(0)))


def read_orth(path):
    """projection -> (ref_gene, q_gene, class)."""
    out = {}
    if not path:
        return out
    with T3C.opener(path) as fh:
        next(fh, None)
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) >= 5 and t[3] != "None":
                out[t[3]] = (t[0], t[2], t[4])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", action="append", required=True,
                    help="NAME=PATTERN with {q} for the query code; order = preference when runs overlap")
    ap.add_argument("--ref-gff", action="append", default=[], help="NAME=GFF of that run's reference annotation")
    ap.add_argument("--name-map", help="TSV ref_gene<TAB>symbol[<TAB>description], all runs")
    ap.add_argument("--include-gcgl", action="store_true", help="also collect genes named glucagon-like (outgroup)")
    ap.add_argument("--status", default="FI,I", help="statuses kept for ortholog projections")
    ap.add_argument("--paralog-status", default="any", help="statuses kept for paralogous projections, or 'any'")
    ap.add_argument("--min-cds", type=int, default=150)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    runs = []
    for r in a.run:
        name, _, pat = r.partition("=")
        if "{q}" not in pat:
            sys.exit(f"[error] --run {r}: pattern needs {{q}}")
        runs.append((name, pat))
    gffs = defaultdict(list)
    for g in a.ref_gff:
        name, _, path = g.partition("=")
        gffs[name].append(path)
    ok_status = set(a.status.split(","))
    ok_para = None if a.paralog_status == "any" else set(a.paralog_status.split(","))
    os.makedirs(a.outdir, exist_ok=True)

    loci = []                                       # one per (run, query, locus)
    for run, pat in runs:
        names = T3C.read_ref_names(gffs.get(run, []), a.name_map)
        qs = queries_of(pat)
        N.info(f"run {run}: {len(qs)} queries, {len(gffs.get(run, []))} reference GFF(s)")
        verdict = {}

        def group_of(tx, gene):
            if gene not in verdict:
                (sym, desc), _ = T3C.name_of(tx, gene, names)
                g = N.name_group(" | ".join([sym, desc]))
                verdict[gene] = g if (g.startswith("exendin") or (a.include_gcgl and g.startswith("GCGL"))) else None
            return verdict[gene]

        n_proj = 0
        for q in qs:
            d = pat.replace("{q}", q)
            bed = T3C.first(d, "query_annotation.bed")
            if not bed:
                N.warn(f"{run}/{q}: no query_annotation.bed")
                continue
            projs = []
            with T3C.opener(bed) as fh:
                for line in fh:
                    t = line.rstrip("\n").split("\t")
                    if len(t) < 6 or line.startswith(("#", "track")):
                        continue
                    tx, gene, flag = T3C.parse_projection(t[3])
                    grp = group_of(tx, gene)
                    if grp:
                        projs.append(dict(projection=t[3], tx=tx, ref_gene=gene, flag=flag, group=grp,
                                          contig=t[0], start=int(t[1]) + 1, end=int(t[2]), strand=t[5]))
            if not projs:
                continue
            status = T3C.read_status(d)
            orth = read_orth(T3C.first(d, "orthology_classification.tsv"))
            fa = T3C.first(d, "nucleotide.fa")
            seqs = T3B.read_fasta(fa, {p["projection"] for p in projs}, label=f"{run}/{q}") if fa else {}
            for p in projs:
                o = orth.get(p["projection"])
                p["status"] = status.get(p["projection"], "NA")
                p["kind"] = ("retro" if p["flag"] in ("retro", "processed_pseudogene") else
                             "paralog" if p["flag"] == "paralog" or (orth and not o) else "ortholog")
                p["orthology_class"] = o[2] if o else ("paralog" if p["kind"] == "paralog" else "NA")
                nt = seqs.get(p["projection"], "")
                p["nt"] = nt[:-3] if nt[-3:] in ("TAA", "TAG", "TGA") else nt
            n_proj += len(projs)
            # cluster this run's projections into loci
            projs.sort(key=lambda p: (p["contig"], p["strand"], p["start"]))
            cur, end = [], -1
            for p in projs + [None]:
                if p is None or (cur and (p["contig"], p["strand"]) != (cur[0]["contig"], cur[0]["strand"])) \
                        or (cur and p["start"] > end):
                    if cur:
                        best = sorted(cur, key=lambda x: (not x["nt"], x["kind"] != "ortholog",
                                                          RANK.get(x["status"], 99), T3B.isoform(x["tx"]),
                                                          -len(x["nt"])))[0]
                        loci.append(dict(best, run=run, query=q, n_proj=len(cur),
                                         lstart=min(x["start"] for x in cur), lend=max(x["end"] for x in cur),
                                         groups=sorted({x["group"] for x in cur}),
                                         ref_genes=sorted({x["ref_genes"] if "ref_genes" in x else x["ref_gene"]
                                                           for x in cur})))
                    cur, end = [], -1
                if p is not None:
                    cur.append(p); end = max(end, p["end"])
        exg = sorted(g for g, v in verdict.items() if v)
        N.info(f"run {run}: {len(exg)} reference gene(s) named exendin{'/GCGL' if a.include_gcgl else ''}: "
               f"{', '.join(exg[:12]) or 'NONE'}; {n_proj} projections")
        if not exg:
            N.warn(f"run {run}: no reference gene named exendin. Its reference may not annotate one, or "
                   f"the names are not resolved (give --ref-gff {run}=...)")

    # ---- the same locus from several runs -> one ----
    order = {r: i for i, (r, _) in enumerate(runs)}
    loci.sort(key=lambda l: (l["query"], l["contig"], l["lstart"]))
    merged = []
    for l in loci:
        hit = next((m for m in merged if m["query"] == l["query"] and m["contig"] == l["contig"]
                    and l["lstart"] <= m["lend"] and l["lend"] >= m["lstart"]), None)
        if hit is None:
            merged.append(dict(l, members=[l]))
            continue
        hit["members"].append(l)
        hit["lstart"], hit["lend"] = min(hit["lstart"], l["lstart"]), max(hit["lend"], l["lend"])
    final = []
    for m in merged:
        best = sorted(m["members"], key=lambda x: (not x["nt"], x["kind"] != "ortholog", RANK.get(x["status"], 99),
                                                  order[x["run"]], -len(x["nt"])))[0]
        final.append(dict(best, lstart=m["lstart"], lend=m["lend"],
                          also_in=",".join(f'{x["run"]}:{x["status"]}' for x in m["members"] if x is not best),
                          groups=sorted({g for x in m["members"] for g in x["groups"]})))

    # ---- filters, names ----
    used, rows = set(), []
    for l in sorted(final, key=lambda l: (l["query"], l["kind"] != "ortholog", l["contig"], l["lstart"])):
        nt = l["nt"]
        aa = T3B.translate(nt)
        why = []
        if l["kind"] == "retro":
            why.append("retrogene")
        elif l["kind"] == "ortholog" and l["status"] not in ok_status:
            why.append(f"status_{l['status']}")
        elif l["kind"] == "paralog" and ok_para is not None and l["status"] not in ok_para:
            why.append(f"paralog_status_{l['status']}")
        if not nt:
            why.append("no_sequence")
        elif len(nt) % 3:
            why.append("cds_not_mult3")
        if nt and len(nt) < a.min_cds:
            why.append(f"cds<{a.min_cds}")
        lab = "_".join(N.group_abbrev(g) for g in l["groups"])
        base = f"{l['query']}_{lab}" + ("p" if l["kind"] == "paralog" else "")
        name, k = base, 0
        while name in used:
            name = base + "bcdefghijklmnopqrstuvwxyz"[k]; k += 1
        used.add(name)
        rows.append(dict(tip=name, query=l["query"], run=l["run"], also_in=l["also_in"] or "-", kind=l["kind"],
                         orthology_class=l["orthology_class"], status=l["status"],
                         group="/".join(l["groups"]), ref_gene=",".join(l.get("ref_genes") or [l["ref_gene"]]), projection=l["projection"],
                         contig=l["contig"], start=l["lstart"], end=l["lend"], strand=l["strand"],
                         cds_bp=len(nt), aa_len=len(aa), internal_stops=aa.count("*"),
                         included="no" if why else "yes", reason=";".join(why) or "-", protein=aa, _nt=nt))

    cols = ["tip", "query", "run", "also_in", "kind", "orthology_class", "status", "group", "ref_gene",
            "projection", "contig", "start", "end", "strand", "cds_bp", "aa_len", "internal_stops",
            "included", "reason", "protein"]
    with open(os.path.join(a.outdir, "candidates.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    keep = [r for r in rows if r["included"] == "yes"]
    with open(os.path.join(a.outdir, "exendin_cds.fa"), "w") as fn, \
         open(os.path.join(a.outdir, "exendin_aa.fa"), "w") as fa, \
         open(os.path.join(a.outdir, "curation.tsv"), "w") as fc:
        fc.write("# Edit 'action' (keep | drop | rename | trim) and 'value', then run 10b_apply_curation.py.\n"
                 "# rename: value = new tip name.  trim: value = first-last amino acid to KEEP (1-based, in the\n"
                 "# unaligned protein), e.g. 43-129 to remove a 42-residue N-terminal extension.\n"
                 "# Write why in 'note'; it is carried into the curation log.\n")
        fc.write("tip\taction\tvalue\tnote\n")
        for r in keep:
            fn.write(f">{r['tip']}\n{r['_nt']}\n")
            fa.write(f">{r['tip']}\n{r['protein'].replace('*', 'X')}\n")
            auto = [r["kind"], r["status"], f"{r['aa_len']} aa"]
            if r["internal_stops"]:
                auto.append(f"{r['internal_stops']} internal stop(s)")
            if not 80 <= r["aa_len"] <= 120:
                auto.append("length outside 80-120 aa: check the model")
            if r["also_in"] != "-":
                auto.append(f"also in {r['also_in']}")
            fc.write(f"{r['tip']}\tkeep\t\t{'; '.join(auto)}\n")

    N.info(f"{len(loci)} exendin loci over all runs -> {len(rows)} after merging runs -> {len(keep)} kept "
           f"from {len({r['query'] for r in keep})} genomes")
    for why, n in sorted(Counter(w for r in rows for w in r["reason"].split(";") if w != "-").items()):
        N.info(f"  left out: {why}: {n}")
    for q in sorted({r["query"] for r in rows}):
        k = [r for r in keep if r["query"] == q]
        N.info(f"  {q}: {len(k)} kept ({sum(r['kind'] == 'paralog' for r in k)} paralogous): "
               + ", ".join(f"{r['tip']} ({r['status']}, {r['aa_len']} aa)" for r in k))


if __name__ == "__main__":
    main()
