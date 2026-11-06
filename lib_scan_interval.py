#!/usr/bin/env python3
"""
Scan a genomic interval for every exendin-family copy, annotated or not.

The tandem array is systematically under-annotated: in HelSus four of six copies
had no gene model, and in HemCap three loci were missing entirely. So presence /
absence has to come from the sequence, not from the GFF.

This runs tblastn with full precursors as queries and classifies every hit by
WHERE IN THE QUERY it falls, which is what tells a signal-peptide exon from a
hormone-domain exon. A hit in the first ~45 residues is the signal/pro region; a
hit after that is the hormone domain. Hits are then clustered into loci, so a
locus with both parts is a candidate gene and a locus with only one is a remnant.

Premature stops in the aligned subject translation are reported (tblastn returns
them as '*' in sseq), as is the first residue of each hormone-domain hit: every
peptide in this family starts with His or Tyr, so anything else is a disabled copy.

Outputs (--out-prefix)
  _hits.tsv   one row per BLAST hit, with part, stops and the translated subject
  _loci.tsv   one row per clustered locus: what parts were found and its state

Usage
  python3 2f_scan_interval.py --genome $GENOME --contig NC_085842.1 \
      --region 155080000-155145000 --query lineage_queries.faa \
      --out-prefix scan_AnoCar

  # re-analyse without re-running BLAST
  python3 2f_scan_interval.py ... --blast-out scan_AnoCar.blast.tsv
"""
import argparse, os, re, subprocess, sys
from collections import defaultdict

MOTIF = re.compile(r"[HY][ASG][DEQ][AGP][ATVIGSM]F[TSI]")
FIELDS = "qseqid qstart qend sstart send pident length evalue bitscore sframe qseq sseq"


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"[error] {cmd[0]} failed:\n{r.stderr.strip()}")
    return r.stdout


def residue_at(qseq, sseq, qstart, target):
    """Subject residue aligned to query position `target`, or '' if not covered."""
    qpos = qstart - 1
    for qc, sc in zip(qseq, sseq):
        if qc != "-":
            qpos += 1
            if qpos == target:
                return sc
    return ""


def compatible(h, x, tol):
    """Can these two hits be exons of one gene? Hits from different queries always can.

    `tol` allows the small query overlap BLAST reports at an exon boundary, where the
    codon is split across the junction and both HSPs claim the shared residue.
    """
    if h["q"] != x["q"]:
        return True
    # hits arrive sorted by genomic position, so h is downstream of x
    return h["qs"] > x["qe"] - tol if h["strand"] == "+" else h["qe"] < x["qs"] + tol


def read_fasta(path):
    out, name, buf = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if name:
        out.append((name, "".join(buf)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", help="indexed genome FASTA (not needed with --blast-out)")
    ap.add_argument("--contig")
    ap.add_argument("--region", help="start-end on the contig")
    ap.add_argument("--query", required=True, help="full precursor proteins, one per lineage")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--blast-out", help="existing tblastn output in the expected -outfmt 6")
    ap.add_argument("--evalue", default="1e-3")
    ap.add_argument("--min-identity", type=float, default=35.0)
    ap.add_argument("--min-length", type=int, default=15)
    ap.add_argument("--cluster-kb", type=float, default=10.0,
                    help="maximum gap between hits of one gene; hits whose QUERY ranges "
                         "overlap are always split, so tandem copies stay separate")
    ap.add_argument("--query-overlap", type=int, default=5,
                    help="query residues two HSPs of one gene may share at an exon boundary")
    ap.add_argument("--default-domain-start", type=int, default=45,
                    help="used for a query with no recognisable hormone-domain motif")
    a = ap.parse_args()

    # where the hormone domain starts in each query
    dstart = {}
    for n, s in read_fasta(a.query):
        m = MOTIF.search(s.upper())
        dstart[n] = (m.start() + 1) if m else a.default_domain_start
    print(f"[info] query domain starts: " +
          ", ".join(f"{n}@{p}" for n, p in dstart.items()), file=sys.stderr)

    offset = 0
    if a.blast_out:
        raw = open(a.blast_out).read()
        if a.region:
            offset = int(a.region.split("-")[0]) - 1
    else:
        if not (a.genome and a.contig):
            sys.exit("--genome and --contig are required unless --blast-out is given")
        region = f"{a.contig}:{a.region}" if a.region else a.contig
        offset = int(a.region.split("-")[0]) - 1 if a.region else 0
        sub = a.out_prefix + ".subject.fa"
        with open(sub, "w") as fh:
            fh.write(run(["samtools", "faidx", a.genome, region]))
        raw = run(["tblastn", "-query", a.query, "-subject", sub,
                   "-outfmt", f"6 {FIELDS}", "-evalue", a.evalue, "-seg", "no"])
        with open(a.out_prefix + ".blast.tsv", "w") as fh:
            fh.write(raw)

    hits = []
    for line in raw.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        q, qs, qe = f[0], int(f[1]), int(f[2])
        ss, se = int(f[3]) + offset, int(f[4]) + offset
        pid, ln = float(f[5]), int(f[6])
        if pid < a.min_identity or ln < a.min_length:
            continue
        d = dstart.get(q, a.default_domain_start)
        part = "N_term" if qe < d else ("hormone_domain" if qs >= d else "spans")
        qseq, sseq = f[10], f[11]
        # every peptide in this family starts His or Tyr: check the residue that
        # aligns to the query's own domain start, for any hit covering it
        first = residue_at(qseq, sseq, qs, d) if qs <= d <= qe else ""
        hits.append(dict(q=q, qs=qs, qe=qe, lo=min(ss, se), hi=max(ss, se),
                         strand="+" if se > ss else "-", pid=pid, ln=ln,
                         ev=f[7], bits=float(f[8]), frame=f[9], sseq=sseq,
                         part=part, stops=sseq.count("*"), first=first))
    if not hits:
        print("[info] no hits passed the filters", file=sys.stderr)

    # Group hits into loci. Several different queries matching the same place are one
    # locus. But two hits FROM THE SAME QUERY belong to the same gene only if they are
    # collinear - a later genomic position must match a later part of that query, as
    # exons of one precursor do. Two tandem copies instead match the same part of the
    # query, which is how they are told apart. Distance alone cannot do it: a gene can
    # span 6 kb while the next copy starts 3 kb later.
    gap = a.cluster_kb * 1000
    hits.sort(key=lambda h: h["lo"])
    loci = []
    for h in hits:
        for L in reversed(loci):
            if h["lo"] - max(x["hi"] for x in L) > gap:
                break
            if L[0]["strand"] != h["strand"]:
                continue
            if not all(compatible(h, x, a.query_overlap) for x in L):
                continue                      # same query, same part -> another copy
            L.append(h)
            break
        else:
            loci.append([h])
            continue
        if not any(h is x for L in loci for x in L):
            loci.append([h])

    with open(a.out_prefix + "_hits.tsv", "w") as fh:
        fh.write("locus\tquery\tq_start\tq_end\tg_start\tg_end\tstrand\tpart\tpident\t"
                 "length\tstops\tfirst_aa\tevalue\tsubject_translation\n")
        for i, L in enumerate(loci, 1):
            for h in sorted(L, key=lambda x: -x["bits"]):
                fh.write(f"{i}\t{h['q']}\t{h['qs']}\t{h['qe']}\t{h['lo']}\t{h['hi']}\t{h['strand']}\t"
                         f"{h['part']}\t{h['pid']:.1f}\t{h['ln']}\t{h['stops']}\t{h['first']}\t"
                         f"{h['ev']}\t{h['sseq']}\n")

    # Two loci that physically overlap in the genome cannot be separate genes: merge
    # them. This catches cases the per-query rule splits because BLAST reported several
    # overlapping HSPs from one query across the same gene.
    changed = True
    while changed:
        changed, merged = False, []
        loci.sort(key=lambda L: min(h["lo"] for h in L))
        for L in loci:
            if merged and min(h["lo"] for h in L) <= max(h["hi"] for h in merged[-1]):
                merged[-1].extend(L); changed = True
            else:
                merged.append(L)
        loci = merged

    rows = []
    for i, L in enumerate(loci, 1):
        parts = {h["part"] for h in L}
        best = max(L, key=lambda h: h["bits"])
        # count stops only in the best hit: a poorly matching secondary query can align
        # through a stop codon in the wrong frame and mark an intact gene as a pseudogene
        stops = best["stops"]
        stops_any = sum(h["stops"] for h in L)
        # judge the domain start from the single best hit that covers it, not from every
        # hit: a weak secondary query can align its domain start anywhere
        dom_hits = [h for h in L if h["first"]]
        bad_start, start_q = [], ""
        if dom_hits:
            b = max(dom_hits, key=lambda h: h["bits"])
            start_q = f"{b['first']}@{b['q']}"
            if b["first"] not in "HY":
                bad_start = [b["first"]]
        if {"N_term", "hormone_domain"} <= parts or "spans" in parts:
            state = "candidate_gene"
        elif parts == {"hormone_domain"}:
            state = "domain_only"
        else:
            state = "N_term_only"
        if stops:
            state += ";premature_stop"
        if bad_start:
            state += f";domain_starts_{bad_start[0]}"
        rows.append([i, min(h["lo"] for h in L), max(h["hi"] for h in L), best["strand"],
                     len(L), "+".join(sorted(parts)), state, best["q"],
                     f"{best['pid']:.1f}", stops, stops_any, start_q])

    with open(a.out_prefix + "_loci.tsv", "w") as fh:
        fh.write("locus\tg_start\tg_end\tstrand\tn_hits\tparts\tstate\tbest_query\t"
                 "best_pident\tstops\tstops_any_hit\tdomain_start\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    print(f"[info] {len(hits)} hits in {len(loci)} loci", file=sys.stderr)
    print("locus\tg_start\tg_end\tstrand\tn_hits\tparts\tstate\tbest_query\tbest_pident\t"
          "stops\tstops_any\tdomain_start", file=sys.stderr)
    for r in rows:
        print("\t".join(map(str, r)), file=sys.stderr)


if __name__ == "__main__":
    main()
