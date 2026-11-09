#!/usr/bin/env python3
"""
Is an exendin-family gene intact, or a pseudogene? Each locus named on the command line
is compared, exon by exon, with the helodermin (exendin-2) precursor (UniProt P04204,
4_reference_peptides.fa EXENDIN2_precursor). Written for the HelSus_B / HelCha_C pair:
the Heloderma copies that match no published peptide (step 5) and sit next to the
exendin-2 gene, between GCGL and EX2.

Per locus
  1. The precursor is searched (sensitive tblastn) against the locus +/- --pad-bp. The
     hits are chained into exons through the best hit at the locus, splice-aware: two
     blocks joined by an intron score better if a GT..AG keeps the frame between them
     and if the intron is short, so a chance hit or the exon of a neighbouring copy (the
     Pro-rich last exons of exendin genes all look alike) does not replace the gene's own
     exon. Precursor stretches no block covers are searched again on their own between
     the neighbouring blocks (gap fill, as 4e) and allowed into the chain. One codon per
     precursor residue (6a), ends extended over mismatches.
  2. Exons are read on the gene's strand. Between two aligned blocks:
       - a gap < --min-intron-bp is inside one exon: genomic length minus 3 x the
         precursor residues it should hold is an indel; not a multiple of 3 = FRAMESHIFT;
       - a longer gap is an intron: a donor (GT, or GC) and acceptor (AG) pair is looked
         for near the block ends that keeps the reading frame (allowing the stretch
         between them to be up to --splice-indel-codons codons longer or shorter than
         the precursor's, reported as an in-frame indel); none = NO CANONICAL SPLICE SITE.
  3. A CDS is assembled from the exons (first exon extended back to precursor residue 1,
     last exon on to its stop codon) and translated. Reported: start codon (ATG?),
     premature stops (also stops inside aligned blocks), frameshifts, CDS length mod 3,
     the protein's identity to the precursor overall and in the signal peptide,
     propeptide and mature peptide, the KR processing site before the mature peptide.
  4. Verdict: 'intact' (ATG, no premature stop, no frameshift, every intron canonical,
     >= 90% of the precursor aligned, a stop where the precursor ends), 'disrupted'
     (any premature stop, frameshift or missing splice site, listed in 'reasons'), else
     'partial' (part of the precursor not found, no disruption).
  5. --pair A,B: Nei-Gojobori dN, dS and dN/dS between two loci, over the codons both
     place against the same precursor residue (so frameshifts do not misalign them; stop
     codons left out; Jukes-Cantor corrected).
     A pseudogene evolves neutrally (dN/dS near 1); give the intact pair as a control.
--toga2-run adds TOGA2's status of the projection at each locus (FI/I/PI/UL/L/M/PG).

Usage
  python3 7a_pseudogene_check.py --species-table ../species.tsv --genomedir ${GENOMEDIR} \
      --curated ${WORK}/reciprocal_blast/reciprocal/curated_loci.tsv \
      --names ${WORK}/nomenclature/names_names.tsv \
      --locus HelSus_B --locus HelCha_C --control HelSus_EX2 --control HelCha_EX1 \
      --pair HelSus_B,HelCha_C --pair HelSus_EX2,HelCha_EX1 --pair HelSus_B,HelSus_EX2 \
      --outdir check

Writes to --outdir
  summary.tsv   one row per locus: verdict, reasons, exons, introns, splice sites,
                frameshifts, stops, identities, TOGA2 status, the aligned proteins
  exons.tsv     the exons: precursor residues, genomic coordinates, splice dinucleotides
  events.tsv    every disruption: locus, kind (stop / frameshift / indel / splice / start),
                precursor position, detail
  residues.tsv  one row per precursor residue and locus (as 6a)
  cds.fna, protein.faa   the assembled CDSs and their translations
  pairs.tsv     dN, dS, dN/dS of each --pair
  reference.tsv the precursor's regions (signal, propeptide, mature)
"""
import argparse, csv, importlib.util, math, os, sys, tempfile
from collections import defaultdict
from itertools import permutations

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402
import lib_gff as G            # noqa: E402

_spec = importlib.util.spec_from_file_location("map6a", os.path.join(HERE, "6a_map_published.py"))
M6 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(M6)
CODE, COMP, translate = M6.CODE, M6.COMP, M6.translate


def read_tsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def nw(a, b, gap=-2):
    """Global alignment, +2 identity / -1 mismatch / gap -2 (end gaps free). -> (aa, ab)."""
    n, m = len(a), len(b)
    S = [[0] * (m + 1) for _ in range(n + 1)]
    T = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        T[i][0] = 1
    for j in range(1, m + 1):
        T[0][j] = 2
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            g1 = 0 if j == m else gap
            g2 = 0 if i == n else gap
            opts = (S[i - 1][j - 1] + (2 if a[i - 1] == b[j - 1] else -1), S[i - 1][j] + g1, S[i][j - 1] + g2)
            k = max(range(3), key=lambda x: opts[x])
            S[i][j], T[i][j] = opts[k], k
    i, j, xa, xb = n, m, [], []
    while i or j:
        k = T[i][j]
        if k == 0:
            xa.append(a[i - 1]); xb.append(b[j - 1]); i -= 1; j -= 1
        elif k == 1:
            xa.append(a[i - 1]); xb.append("-"); i -= 1
        else:
            xa.append("-"); xb.append(b[j - 1]); j -= 1
    return "".join(xa[::-1]), "".join(xb[::-1])


# ---------------- Nei-Gojobori ----------------
def _syn_sites(c):
    if CODE.get(c, "X") in "*X":
        return None
    s = 0.0
    for p in range(3):
        for b in "ACGT":
            if b == c[p]:
                continue
            m = c[:p] + b + c[p + 1:]
            if CODE[m] == "*":
                continue
            s += (CODE[m] == CODE[c]) / 3.0
    return s


def _diffs(c1, c2):
    """(synonymous, nonsynonymous) differences averaged over pathways without stops."""
    pos = [i for i in range(3) if c1[i] != c2[i]]
    if not pos:
        return 0.0, 0.0
    tot_s = tot_n = 0.0
    npath = 0
    for order in permutations(pos):
        cur, s, n, ok = c1, 0, 0, True
        for p in order:
            nxt = cur[:p] + c2[p] + cur[p + 1:]
            if CODE[nxt] == "*":
                ok = False; break
            if CODE[nxt] == CODE[cur]:
                s += 1
            else:
                n += 1
            cur = nxt
        if ok:
            tot_s += s; tot_n += n; npath += 1
    if not npath:
        return 0.0, float(len(pos))
    return tot_s / npath, tot_n / npath


def ng86(cds1, cds2, aln1, aln2):
    """dN/dS of two CDSs whose proteins are aligned as aln1/aln2."""
    i = j = 0
    S = Nn = Sd = Nd = 0.0
    used = 0
    for x, y in zip(aln1, aln2):
        c1 = cds1[3 * i:3 * i + 3] if x != "-" else None
        c2 = cds2[3 * j:3 * j + 3] if y != "-" else None
        if x != "-":
            i += 1
        if y != "-":
            j += 1
        if not c1 or not c2 or len(c1) < 3 or len(c2) < 3 or "N" in c1 + c2:
            continue
        s1, s2 = _syn_sites(c1), _syn_sites(c2)
        if s1 is None or s2 is None:
            continue
        s = (s1 + s2) / 2
        S += s; Nn += 3 - s
        sd, nd = _diffs(c1, c2)
        Sd += sd; Nd += nd; used += 1
    if not used or S == 0 or Nn == 0:
        return dict(codons=used)
    pS, pN = Sd / S, Nd / Nn

    def jc(p):
        x = 1 - 4 * p / 3
        return -0.75 * math.log(x) if x > 0 else float("inf")
    dS, dN = jc(pS), jc(pN)
    w = dN / dS if dS not in (0, float("inf")) else float("nan")
    return dict(codons=used, S=S, N=Nn, Sd=Sd, Nd=Nd, pS=pS, pN=pN, dS=dS, dN=dN, w=w)


# ---------------- splice-aware chaining ----------------
def splice_pair(O, L, o1, o0, qg, k_max):
    """Donor/acceptor joining an exon block ending at oriented o1 to one starting at o0,
    with qg precursor residues unaligned between them. The intron must keep the frame:
    the coding bases left between the blocks (donor side + acceptor side) = 3 x (qg + k),
    k = -k_max..k_max codons (a diverged stretch may be a few codons longer or shorter).
    -> (score, first intron base, last intron base, donor, acceptor, k) or None."""
    best = None
    q = max(qg, 0)
    for k in range(-k_max, k_max + 1):
        tot = 3 * (qg + k)
        if tot < -36:
            continue
        for do in range(-18, 3 * min(q + max(k, 0), 10) + 19):
            ao = tot - do
            if ao < -18 or ao > 3 * min(q + max(k, 0), 10) + 18:
                continue
            d = o1 + 1 + do
            ac = o0 - 1 - ao
            if d < 0 or ac >= L or ac - d < 20:
                continue
            dn, an = O[d:d + 2], O[ac - 1:ac + 1]
            if an != "AG" or dn not in ("GT", "GC"):
                continue
            sc = (dn != "GT") * 100 + 3 * abs(k) + (abs(do) + abs(ao)) / 10
            if best is None or sc < best[0]:
                best = (sc, d, ac, dn, an, k)
    return best


def dp_chain(cands, anchor, strand, O, L, ori, a):
    """Chain of HSP blocks through the anchor, scored as  sum(bits) - junction penalties:
      inside an exon (gap < --min-intron-bp): 0 for no indel, 2 for an in-frame indel,
        --frameshift-pen for a frameshift;
      intron: --intron-pen + 1 bit per kb of intron + 3 per codon of length difference,
        + --nosplice-pen if no GT/GC..AG keeps the frame.
    So the gene's own exons, joined by real introns, beat chance hits and the exons of a
    neighbouring copy (e.g. the Pro-rich last exon of another exendin gene)."""
    hs = sorted(cands, key=lambda h: (h["qs"], h["qe"]))
    for h in hs:
        h["o0"], h["o1"] = sorted([ori(h["lo"]), ori(h["hi"])])
    n = len(hs)
    cache = {}

    def pen(i, j):
        if (i, j) in cache:
            return cache[(i, j)]
        A, B = hs[i], hs[j]
        out = None
        if A["qs"] < B["qs"] and A["qe"] < B["qe"] and A["qe"] - B["qs"] + 1 <= 8:
            g = B["o0"] - A["o1"] - 1
            qg = B["qs"] - A["qe"] - 1
            if -30 <= g < a.min_intron_bp:
                indel = g - 3 * qg
                out = 0.0 if indel == 0 else (2.0 if indel % 3 == 0 else a.frameshift_pen)
            elif a.min_intron_bp <= g <= a.max_intron_bp:
                sp = splice_pair(O, L, A["o1"], B["o0"], qg, a.splice_indel_codons)
                out = a.intron_pen + g / 1000.0 + (sp[0] if sp and sp[0] < 100 else
                                                   (3 + sp[0] - 100 if sp else a.nosplice_pen))
        cache[(i, j)] = out
        return out
    F, Fp = [h["bits"] for h in hs], [None] * n
    for j in range(n):
        for i in range(j):
            p_ = pen(i, j)
            if p_ is not None and F[i] - p_ + hs[j]["bits"] > F[j]:
                F[j], Fp[j] = F[i] - p_ + hs[j]["bits"], i
    B_, Bp = [h["bits"] for h in hs], [None] * n
    for i in range(n - 1, -1, -1):
        for j in range(i + 1, n):
            p_ = pen(i, j)
            if p_ is not None and B_[j] - p_ + hs[i]["bits"] > B_[i]:
                B_[i], Bp[i] = B_[j] - p_ + hs[i]["bits"], j
    ai = next(i for i, h in enumerate(hs) if h is anchor)
    left, i = [], ai
    while i is not None:
        left.append(hs[i]); i = Fp[i]
    right, i = [], Bp[ai]
    while i is not None:
        right.append(hs[i]); i = Bp[i]
    return left[::-1] + right


# ---------------- per-locus analysis ----------------
def fill_gaps(ch, prec, plus, w0, w1, contig, strand, a, tmp):
    """Precursor stretches (>= 5 residues) no block covers are searched again on their own,
    in the stretch of genome where they must lie: before the first block, between two
    blocks, or after the last (as 4e's gap fill). A short diverged exon scores too little
    inside a full-length search to join the chain; alone, in its own window, it is found.
    Accepted: best hit on the strand with >= --fill-bits, >= 50% identity, covering
    >= min(8, stretch) residues."""
    ch = sorted(ch, key=lambda h: h["qs"])
    if not ch:
        return ch
    gaps = []
    if ch[0]["qs"] - 1 >= 5:
        gaps.append((1, ch[0]["qs"] - 1, None, ch[0]))
    for x, y in zip(ch, ch[1:]):
        if y["qs"] - x["qe"] - 1 >= 5:
            gaps.append((x["qe"] + 1, y["qs"] - 1, x, y))
    if len(prec) - ch[-1]["qe"] >= 5:
        gaps.append((ch[-1]["qe"] + 1, len(prec), ch[-1], None))
    added = []
    for q0, q1, left, right in gaps:
        # genomic interval between the neighbouring blocks, on the gene's strand
        if strand == "+":
            lo = left["hi"] + 1 if left else right["lo"] - a.max_intron_bp
            hi = right["lo"] - 1 if right else left["hi"] + a.max_intron_bp
        else:
            lo = right["hi"] + 1 if right else left["lo"] - a.max_intron_bp
            hi = left["lo"] - 1 if left else right["hi"] + a.max_intron_bp
        lo, hi = max(lo, w0), min(hi, w1)
        if hi - lo < 3 * (q1 - q0 + 1) // 2:
            continue
        sub = plus[lo - w0:hi - w0 + 1]
        sf, pf = os.path.join(tmp, "g.fa"), os.path.join(tmp, "g.faa")
        open(sf, "w").write(f">g\n{sub}\n")
        open(pf, "w").write(f">s\n{prec[q0 - 1:q1]}\n")
        hs = M6.parse(M6.run(["tblastn", "-query", pf, "-subject", sf, "-seg", "no", "-comp_based_stats", "0",
                              "-word_size", "2", "-evalue", "1000", "-outfmt", M6.OUTFMT]), offset=(contig, lo))
        need = min(8, q1 - q0 + 1)
        hs = [h for h in hs if h["strand"] == strand and h["bits"] >= a.fill_bits and h["pident"] >= 50
              and h["qe"] - h["qs"] + 1 >= need]
        if not hs:
            continue
        h = max(hs, key=lambda h: h["bits"])
        h["qs"] += q0 - 1; h["qe"] += q0 - 1; h["fill"] = True
        added.append(h)
    return sorted(ch + added, key=lambda h: h["qs"])



def analyse(name, rec, files, prec, regions, a, tmp):
    contig = rec["contig"]
    clen = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(files["fasta"] + ".fai")}
    w0 = max(1, int(rec["g_start"]) - a.pad_bp)
    w1 = min(clen[contig], int(rec["g_end"]) + a.pad_bp)
    plus = M6.faidx(files["fasta"], contig, w0, w1)
    sf, pf = os.path.join(tmp, "w.fa"), os.path.join(tmp, "p.faa")
    open(sf, "w").write(f">w\n{plus}\n")
    open(pf, "w").write(f">prec\n{prec}\n")
    hits = M6.parse(M6.run(["tblastn", "-query", pf, "-subject", sf, "-seg", "no", "-comp_based_stats", "0",
                            "-word_size", "2", "-evalue", "1000", "-outfmt", M6.OUTFMT]), offset=(contig, w0))
    strand = rec.get("strand") if rec.get("strand") in ("+", "-") else None
    lo_c, hi_c = int(rec["g_start"]), int(rec["g_end"])
    near = [h for h in hits if h["hi"] >= lo_c - 2000 and h["lo"] <= hi_c + 2000 and (not strand or h["strand"] == strand)]
    if not near:
        return dict(name=name, verdict="not found", reasons="precursor not found at the locus"), [], [], [], "", ""
    anchor = max(near, key=lambda h: h["bits"])
    strand = anchor["strand"]
    L = w1 - w0 + 1
    O = plus if strand == "+" else plus.translate(COMP)[::-1]

    def ori(p):                                         # genomic -> oriented 0-based
        return p - w0 if strand == "+" else w1 - p
    cands = [h for h in hits if h["strand"] == strand and
             (h is anchor or (h["bits"] >= a.min_exon_bits and h["qe"] - h["qs"] + 1 >= 4))]
    ch = dp_chain(cands, anchor, strand, O, L, ori, a)
    filled = [h for h in fill_gaps(ch, prec, plus, w0, w1, contig, strand, a, tmp) if h.get("fill")]
    if filled:                                          # let the re-searched stretches compete
        ch = dp_chain(cands + filled, anchor, strand, O, L, ori, a)
    rows, _, _ = M6.residue_map(ch, prec)
    rows, _ = M6.extend_ends(rows, files["fasta"], contig, strand, clen[contig], max_ext=6)

    def gen(o):                                         # oriented 0-based -> genomic
        return w0 + o if strand == "+" else w1 - o
    # aligned blocks (one per chain exon), oriented
    segs = defaultdict(list)
    for r in rows:
        if r["codon_start"] != "" and r["status"] != "gap":
            o = ori(int(r["codon_start"])) if strand == "+" else ori(int(r["codon_end"]))
            segs[r["exon"]].append((r["qpos"], o))
    blocks = []
    for e, v in segs.items():
        v.sort()
        blocks.append(dict(q0=v[0][0], q1=v[-1][0], o0=v[0][1], o1=v[-1][1] + 2))
    blocks.sort(key=lambda b: b["q0"])
    events = []
    for r in rows:
        if r["status"] == "stop":
            events.append(dict(kind="stop", qpos=r["qpos"], detail=f"stop codon aligned to precursor {prec[r['qpos'] - 1]}{r['qpos']}"))
    # merge blocks into exons; check junctions
    exons = [dict(q0=blocks[0]["q0"], q1=blocks[0]["q1"], o0=blocks[0]["o0"], o1=blocks[0]["o1"])] if blocks else []
    introns = []
    for b in blocks[1:]:
        p = exons[-1]
        g = b["o0"] - p["o1"] - 1
        qg = b["q0"] - p["q1"] - 1
        if g < a.min_intron_bp:
            indel = g - 3 * qg
            if indel % 3:
                events.append(dict(kind="frameshift", qpos=p["q1"],
                                   detail=f"{indel:+d} bp between precursor {p['q1']} and {b['q0']} (not a multiple of 3)"))
            elif indel:
                events.append(dict(kind="indel", qpos=p["q1"],
                                   detail=f"in-frame {indel:+d} bp ({indel // 3:+d} codons) after precursor {p['q1']}"))
            p["q1"], p["o1"] = b["q1"], b["o1"]
            continue
        # intron: donor/acceptor that keep the frame
        best = splice_pair(O, L, p["o1"], b["o0"], qg, a.splice_indel_codons)
        if best:
            _, d, ac, dn, an, k = best
            if k:
                events.append(dict(kind="indel", qpos=p["q1"],
                                   detail=f"in-frame {k:+d} codon(s) between precursor {p['q1']} and {b['q0']} "
                                          "(diverged stretch of another length)"))
            introns.append(dict(after_q=p["q1"], donor=dn, acceptor=an, canonical="yes", length=ac - d + 1))
            if dn == "GC":
                events.append(dict(kind="splice", qpos=p["q1"], detail="GC donor (non-canonical but functional in ~1% of introns)"))
            p["o1"] = d - 1
            exons.append(dict(q0=b["q0"], q1=b["q1"], o0=ac + 1, o1=b["o1"]))
        else:
            dn, an = O[p["o1"] + 1:p["o1"] + 3], O[b["o0"] - 2:b["o0"]]
            introns.append(dict(after_q=p["q1"], donor=dn, acceptor=an, canonical="no", length=g))
            events.append(dict(kind="splice", qpos=p["q1"],
                               detail=f"no GT/GC..AG pair keeps the frame between precursor {p['q1']} and {b['q0']} "
                                      f"(block ends read {dn}..{an})"))
            exons.append(dict(q0=b["q0"], q1=b["q1"], o0=b["o0"], o1=b["o1"]))
    if not exons:
        return dict(name=name, verdict="not found", reasons="no aligned block"), [], events, rows, "", ""
    # extend the first exon back to residue 1 and the last on to the stop codon, if only a
    # few residues are missing there (more = an exon not found: not extended through the intron)
    head_missing = exons[0]["q0"] - 1
    tail_missing = len(prec) - exons[-1]["q1"]
    if head_missing <= a.max_end_ext:
        exons[0]["o0"] = max(0, exons[0]["o0"] - 3 * head_missing)
    if tail_missing <= a.max_end_ext:
        exons[-1]["o1"] = min(L - 1, exons[-1]["o1"] + 3 * tail_missing + 3)
    cds = "".join(O[e["o0"]:e["o1"] + 1] for e in exons)
    prot = translate(cds)
    start = cds[:3] if head_missing <= a.max_end_ext else "-"
    if start not in ("ATG", "-"):
        events.append(dict(kind="start", qpos=1, detail=f"no ATG at precursor M1 (codon {start})"))
    xa, xb = nw(prec, prot.rstrip("*") if prot.endswith("*") and "*" not in prot[:-1] else prot)
    # where the assembled ORF ends, in precursor coordinates (stops after a frameshift are
    # out of frame and not counted separately)
    first_stop = prot.find("*")
    orf_end_q = None
    if first_stop >= 0 and first_stop < len(prot) - 1:
        qi = ti = 0
        for x, y in zip(xa, xb):
            if x != "-":
                qi += 1
            if y != "-":
                ti += 1
                if ti == first_stop + 1:
                    orf_end_q = qi
                    break
        events.append(dict(kind="orf_end", qpos=orf_end_q or 0,
                           detail=f"assembled ORF ends after {first_stop} codons (at precursor ~{orf_end_q})"))
    if len(cds) % 3 and not any(e["kind"] == "frameshift" for e in events):
        events.append(dict(kind="frameshift", qpos=exons[-1]["q1"], detail=f"CDS length {len(cds)} not a multiple of 3"))
    term_ok = first_stop == len(prot) - 1 if tail_missing <= a.max_end_ext else None
    # identities by region, from the codon-by-codon mapping (independent of the frame of the
    # assembled CDS, so a frameshift does not hide how similar the rest of the gene is)
    pos_id = {r["qpos"]: r["status"] == "identical" for r in rows}
    def ident(r0, r1):
        n = r1 - r0 + 1
        k = sum(1 for q in range(r0, r1 + 1) if pos_id.get(q))
        return k, n
    aligned = sum(1 for r in rows if r["codon_start"] != "")
    cover = aligned / len(prec)
    stops = [e for e in events if e["kind"] == "stop"]
    fs = [e for e in events if e["kind"] == "frameshift"]
    spl = [e for e in events if e["kind"] == "splice" and "no GT" in e["detail"]]
    reasons = []
    if stops:
        reasons.append("premature stop at " + ",".join(str(e["qpos"]) for e in stops))
    if fs:
        reasons.append("frameshift after " + ",".join(str(e["qpos"]) for e in fs))
    if spl:
        reasons.append("no splice site after " + ",".join(str(e["qpos"]) for e in spl))
    if start not in ("ATG", "-"):
        reasons.append(f"start codon {start}")
    if head_missing > a.max_end_ext:
        reasons.append(f"precursor 1-{head_missing} not found (first exon missing or diverged)")
    if tail_missing > a.max_end_ext:
        reasons.append(f"precursor {exons[-1]['q1'] + 1}-{len(prec)} not found (last exon missing or diverged)")
    if term_ok is False and not stops and not fs:
        reasons.append("no stop codon where the precursor ends")
    orf = [e for e in events if e["kind"] == "orf_end"]
    if orf and not (stops or fs or spl):
        reasons.append(orf[0]["detail"])
    if stops or fs or spl or orf:
        verdict = "disrupted"
    elif cover >= 0.9 and start == "ATG" and term_ok and not spl:
        verdict = "intact"
    else:
        verdict = "partial"
        if cover < 0.9:
            reasons.append(f"only {aligned}/{len(prec)} precursor residues found")
    m0, m1 = regions["mature"]
    byq = {r["qpos"]: r for r in rows}                  # codons on precursor m0-2, m0-1 (frame-independent)
    kr = "".join((byq[q]["gres"] if byq.get(q) and byq[q]["codon_start"] != "" else "·") for q in (m0 - 2, m0 - 1))
    sig_k, sig_n = ident(*regions["signal"])
    pro_k, pro_n = ident(*regions["propeptide"])
    mat_k, mat_n = ident(m0, m1)
    all_k, all_n = ident(1, len(prec))
    ex_rows = []
    for k, e in enumerate(exons, 1):
        g0, g1 = sorted([gen(e["o0"]), gen(e["o1"])])
        intr = introns[k - 1] if k - 1 < len(introns) else None
        ex_rows.append(dict(locus=name, exon=k, q_start=e["q0"] if k > 1 else 1,
                            q_end=e["q1"] if k < len(exons) else len(prec), contig=contig, strand=strand,
                            g_start=g0, g_end=g1, length=e["o1"] - e["o0"] + 1,
                            donor_after=intr["donor"] if intr else "", acceptor_after=intr["acceptor"] if intr else "",
                            splice_ok=intr["canonical"] if intr else "",
                            intron_after=intr["length"] if intr else ""))
    summ = dict(name=name, species=rec["species"], step4_name=rec["name"], contig=contig, strand=strand,
                start=min(r["g_start"] for r in ex_rows), end=max(r["g_end"] for r in ex_rows),
                verdict=verdict, reasons="; ".join(reasons) or "-", exons=len(exons),
                introns=",".join(f"{i['length']}" for i in introns) or "-",
                splice=",".join(f"{i['donor']}..{i['acceptor']}" for i in introns) or "-",
                start_codon=start, cds_len=len(cds), protein_len=len(prot.rstrip("*")),
                first_stop=(first_stop + 1) if first_stop >= 0 else "none", aligned=aligned,
                precursor_len=len(prec), identity=f"{all_k}/{all_n}", signal_id=f"{sig_k}/{sig_n}",
                propeptide_id=f"{pro_k}/{pro_n}", mature_id=f"{mat_k}/{mat_n}", kr_site=kr or "-",
                stops=len(stops), frameshifts=len(fs), aln_ref=xa, aln_target=xb)
    codons = {}
    for r in rows:
        if r["codon_start"] != "" and r["status"] not in ("gap", "stop"):
            o = ori(int(r["codon_start"])) if strand == "+" else ori(int(r["codon_end"]))
            codons[r["qpos"]] = O[o:o + 3]
    summ["_codons"] = codons
    return summ, ex_rows, events, rows, cds, prot


def toga_names(fasta, extra):
    """TOGA2 contig name -> genome contig: --contig-map, the NCBI sequence report (chr2 ->
    NC_...), TOGA2's renaming (JANEZZ010002025__1 -> JANEZZ010002025.1)."""
    import json, re
    fai = {l.split("\t")[0] for l in open(fasta + ".fai")}
    m = {}
    rp = os.path.join(os.path.dirname(fasta), "sequence_report.jsonl")
    if os.path.exists(rp):
        for line in open(rp):
            d = json.loads(line)
            acc = d.get("refseqAccession") if d.get("refseqAccession") in fai else d.get("genbankAccession")
            for k in ("chrName", "sequenceName", "genbankAccession", "refseqAccession"):
                if d.get(k) and acc:
                    m[d[k]] = acc
                    if k == "chrName":
                        m["chr" + d[k]] = acc
    for cm in extra:
        t, _, n = cm.partition("=")
        m[t.strip()] = n.strip()

    def f(t):
        if t in m:
            return m[t]
        if t in fai:
            return t
        for c in (t.replace("__", "."), re.sub(r"_(\d+)$", r".\1", t)):
            if c in fai:
                return c
        return t
    return f


def regions_of(prec, mature):
    i = prec.find(mature[:20])
    if i < 0:
        sys.exit("[error] mature peptide not found in the precursor")
    m0 = i + 1
    m1 = min(len(prec), m0 + len(mature) - 1)
    sig = (1, min(22, m0 - 1))
    return {"signal": sig, "propeptide": (sig[1] + 1, m0 - 1), "mature": (m0, m1)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--genomedir", required=True)
    ap.add_argument("--curated", required=True, help="step 4 curated_loci.tsv")
    ap.add_argument("--names", help="step 5 names_names.tsv (lets --locus use the new names)")
    ap.add_argument("--locus", action="append", default=[], help="locus to test (new or step-4 name); repeatable")
    ap.add_argument("--control", action="append", default=[], help="locus expected intact; repeatable")
    ap.add_argument("--pair", action="append", default=[], metavar="A,B", help="dN/dS between two loci; repeatable")
    ap.add_argument("--reference", default=os.path.join(HERE, "4_reference_peptides.fa"))
    ap.add_argument("--precursor-id", default="EXENDIN2_precursor")
    ap.add_argument("--mature-id", default="EXENDIN2_mature")
    ap.add_argument("--signal-len", type=int, default=22, help="signal peptide length of the precursor [22]")
    ap.add_argument("--pad-bp", type=int, default=10000)
    ap.add_argument("--max-intron-bp", type=int, default=10000)
    ap.add_argument("--min-intron-bp", type=int, default=60)
    ap.add_argument("--max-exons", type=int, default=6)
    ap.add_argument("--min-exon-bits", type=float, default=16.0,
                    help="bits an aligned block other than the strongest needs [16]")
    ap.add_argument("--frameshift-pen", type=float, default=15.0, help="chaining penalty for a frameshift [15]")
    ap.add_argument("--intron-pen", type=float, default=4.0, help="chaining penalty per intron [4], plus 1 per kb")
    ap.add_argument("--nosplice-pen", type=float, default=25.0,
                    help="chaining penalty for an intron with no frame-keeping GT/GC..AG [25]")
    ap.add_argument("--splice-indel-codons", type=int, default=4,
                    help="a diverged stretch across an intron may be this many codons longer or shorter [4]")
    ap.add_argument("--fill-bits", type=float, default=12.0,
                    help="bits a re-searched precursor stretch needs to be accepted as an exon [12]")
    ap.add_argument("--max-end-ext", type=int, default=6,
                    help="extend the CDS over at most this many unaligned precursor residues at each end [6]; the last exon is ~10 residues, so a longer stretch means it was not found")
    ap.add_argument("--toga2-run", help="${TOGADIR}/target_HelSus: add TOGA2's status of the projection at each locus")
    ap.add_argument("--contig-map", action="append", default=[], metavar="TOGA2=GENOME",
                    help="TOGA2 contig name -> genome contig, if the automatic translation fails")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    ref = N.read_fasta(a.reference)
    pk = next((k for k in ref if k.split("|")[0] == a.precursor_id), None)
    mk = next((k for k in ref if k.split("|")[0] == a.mature_id), None)
    if not pk or not mk:
        sys.exit(f"[error] {a.precursor_id} / {a.mature_id} not in {a.reference}")
    prec, mature = ref[pk][1].upper(), ref[mk][1].upper()
    regions = regions_of(prec, mature)
    regions["signal"] = (1, min(a.signal_len, regions["mature"][0] - 1))
    regions["propeptide"] = (regions["signal"][1] + 1, regions["mature"][0] - 1)
    N.info(f"reference {pk.split('|')[0]} ({len(prec)} aa): signal {regions['signal']}, propeptide "
           f"{regions['propeptide']}, mature {regions['mature']}")

    table = {r["code"]: r for r in N.read_species_table(a.species_table)}
    cur = {r["name"]: r for r in read_tsv(a.curated)}
    new2old = {}
    if a.names:
        for r in read_tsv(a.names):
            new2old[r["new_name"]] = r["old_name"]
    old2new = {v: k for k, v in new2old.items()}

    def resolve(nm):
        if nm in cur:
            return old2new.get(nm, nm), cur[nm]
        if nm in new2old and new2old[nm] in cur:
            return nm, cur[new2old[nm]]
        sys.exit(f"[error] {nm}: not a curated locus" + (" (or a step-5 name)" if a.names else " (give --names?)"))
    todo = [(resolve(x), "test") for x in a.locus] + [(resolve(x), "control") for x in a.control]
    if not todo:
        sys.exit("[error] give --locus and/or --control")

    toga = {}
    tmp = tempfile.mkdtemp(prefix="7a_")
    S, EX, EV, RES, CDS, PROT, COD = [], [], [], [], {}, {}, {}
    for (name, rec), role in todo:
        code = rec["species"]
        files = N.find_genome_files(a.genomedir, code, table[code]["accession"])
        if not os.path.exists(files["fasta"] + ".fai"):
            M6.run(["samtools", "faidx", files["fasta"]])
        summ, ex, ev, rows, cds, prot = analyse(name, rec, files, prec, regions, a, tmp)
        summ["role"] = role
        if a.toga2_run and summ.get("contig"):
            q = os.path.join(a.toga2_run, "query_" + code)
            if code not in toga:
                toga[code] = G.load_toga2_query(q) if os.path.isdir(q) else {}
            names = toga_names(files["fasta"], a.contig_map)
            st = []
            for tc, gs in toga[code].items():
                if names(tc) != summ["contig"]:
                    continue
                st += [f"{g['name']}:{g['status']}" for g in gs if g["end"] >= int(summ["start"]) and g["start"] <= int(summ["end"])]
            summ["toga2"] = ",".join(st) or ("none" if toga[code] else "no query dir")
        S.append(summ); EX += ex
        EV += [dict(locus=name, **e) for e in ev]
        RES += [dict(locus=name, **r) for r in rows]
        CDS[name], PROT[name] = cds, prot
        COD[name] = summ.pop("_codons", {})
        N.info(f"{name} ({role}): {summ['verdict']} | {summ.get('reasons', '')} | exons {summ.get('exons', '-')}, "
               f"splice {summ.get('splice', '-')}, identity {summ.get('identity', '-')}, mature {summ.get('mature_id', '-')}"
               + (f", TOGA2 {summ['toga2']}" if "toga2" in summ else ""))
    pairs = []
    for pr in a.pair:
        x, y = [z.strip() for z in pr.split(",")]
        x = resolve(x)[0]; y = resolve(y)[0]
        if x not in COD or y not in COD:
            N.warn(f"pair {pr}: both loci must be in --locus/--control"); continue
        qs = sorted(set(COD[x]) & set(COD[y]))
        r = ng86("".join(COD[x][q] for q in qs), "".join(COD[y][q] for q in qs), "A" * len(qs), "A" * len(qs))
        r = {k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()}
        pairs.append(dict(a=x, b=y, **r))
        N.info(f"dN/dS {x} vs {y}: " + (f"dN {r.get('dN')} dS {r.get('dS')} w {r.get('w')} ({r['codons']} codons)"
                                        if "w" in r else f"not computable ({r['codons']} codons)"))

    def write(fn, cols, rows):
        with open(os.path.join(a.outdir, fn), "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    write("summary.tsv", ["name", "role", "species", "step4_name", "contig", "strand", "start", "end", "verdict", "reasons",
                          "exons", "introns", "splice", "start_codon", "cds_len", "protein_len", "first_stop", "aligned",
                          "precursor_len", "identity", "signal_id", "propeptide_id", "mature_id", "kr_site", "stops",
                          "frameshifts", "toga2", "aln_ref", "aln_target"], S)
    write("exons.tsv", ["locus", "exon", "q_start", "q_end", "contig", "strand", "g_start", "g_end", "length",
                        "donor_after", "acceptor_after", "splice_ok", "intron_after"], EX)
    write("events.tsv", ["locus", "kind", "qpos", "detail"], EV)
    write("residues.tsv", ["locus", "qpos", "qres", "exon", "codon_start", "codon_end", "gres", "status", "ins_after",
                           "ext"], RES)
    write("pairs.tsv", ["a", "b", "codons", "S", "N", "Sd", "Nd", "pS", "pN", "dS", "dN", "w"], pairs)
    write("reference.tsv", ["region", "start", "end", "sequence"],
          [dict(region=k, start=v[0], end=v[1], sequence=prec[v[0] - 1:v[1]]) for k, v in regions.items()])
    with open(os.path.join(a.outdir, "cds.fna"), "w") as fh:
        for k, v in CDS.items():
            if v:
                fh.write(f">{k}\n{v}\n")
    with open(os.path.join(a.outdir, "protein.faa"), "w") as fh:
        fh.write(f">{a.precursor_id}\n{prec}\n")
        for k, v in PROT.items():
            if v:
                fh.write(f">{k}\n{v}\n")
    N.info(f"wrote {a.outdir}/summary.tsv, exons.tsv, events.tsv, pairs.tsv")


if __name__ == "__main__":
    main()
