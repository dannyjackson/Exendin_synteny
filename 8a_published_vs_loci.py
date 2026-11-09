#!/usr/bin/env python3
"""
Every published exendin against every exendin-family gene of the Heloderma genomes, to
support which species each peptide was purified from.

The published peptides (8_published_peptides.fa: exendin-1 = helospectin I, exendin-2 =
helodermin, exendin-3, exendin-4) were purified from venom pooled from lizards sold as
H. suspectum or H. horridum. A peptide encoded EXACTLY by a gene of one species, and only
approximately by every gene of the other, was most likely purified from the first. Each
peptide is therefore mapped, codon by codon, onto each curated exendin-family locus
(step 4, named by step 5) of each genome (--species, default HelSus,HelCha), not onto the
best locus of the genome only.

Per peptide and locus: sensitive tblastn of the peptide against the locus +/- --pad-bp,
chained through the best hit inside the locus with 7a's splice-aware chaining and gap
fill (the mature peptides of exendin-3/-4 span two exons: the last one is ~10 residues),
one codon per residue (6a), ends extended over mismatches; the genome residues just
before residue 1 (processing site, KR) and just after the last residue (Gly = amidation
signal) are read from the genome.

Usage
  python3 8a_published_vs_loci.py --species-table ../species.tsv --genomedir ${GENOMEDIR} \
      --curated ${WORK}/reciprocal_blast/reciprocal/curated_loci.tsv \
      --names ${WORK}/nomenclature/names_names.tsv --species HelSus,HelCha --outdir match

Writes to --outdir
  matches.tsv    one row per peptide x locus: identical / length, exons, before, after,
                 coordinates, rank of the locus for that peptide within its genome
                 (best_in_species = yes for the best), and the best genome for the peptide
  best.tsv       per peptide: best locus and identity in each genome, and the margin
  residues.tsv   one row per peptide residue and locus (as 6a)
"""
import argparse, csv, importlib.util, os, sys, tempfile
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402


def _load(name, fn):
    sp = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


M6 = _load("map6a", "6a_map_published.py")
M7 = _load("check7a", "7a_pseudogene_check.py")


def read_tsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def map_peptide(pep, rec, files, clen, plus, w0, w1, a, tmp):
    """-> (rows, identical, exons, strand, before, after) or None if no hit at the locus."""
    contig = rec["contig"]
    sf, pf = os.path.join(tmp, "w.fa"), os.path.join(tmp, "p.faa")
    open(sf, "w").write(f">w\n{plus}\n")
    open(pf, "w").write(f">p\n{pep}\n")
    hits = M6.parse(M6.run(["tblastn", "-query", pf, "-subject", sf, "-seg", "no", "-comp_based_stats", "0",
                            "-word_size", "2", "-evalue", "1000", "-outfmt", M6.OUTFMT]), offset=(contig, w0))
    lo_c, hi_c = int(rec["g_start"]), int(rec["g_end"])
    st = rec.get("strand") if rec.get("strand") in ("+", "-") else None
    inside = [h for h in hits if h["hi"] >= lo_c - a.inside_bp and h["lo"] <= hi_c + a.inside_bp
              and (not st or h["strand"] == st) and h["bits"] >= a.min_anchor_bits]
    if not inside:
        return None
    anchor = max(inside, key=lambda h: h["bits"])
    strand = anchor["strand"]
    L = w1 - w0 + 1
    O = plus if strand == "+" else plus.translate(M6.COMP)[::-1]

    def ori(p):
        return p - w0 if strand == "+" else w1 - p
    cands = [h for h in hits if h["strand"] == strand and
             (h is anchor or (h["bits"] >= a.min_exon_bits and h["qe"] - h["qs"] + 1 >= 4))]
    ch = M7.dp_chain(cands, anchor, strand, O, L, ori, a)
    filled = [h for h in M7.fill_gaps(ch, pep, plus, w0, w1, contig, strand, a, tmp) if h.get("fill")]
    if filled:
        ch = M7.dp_chain(cands + filled, anchor, strand, O, L, ori, a)
    rows, _, _ = M6.residue_map(ch, pep)
    rows, _ = M6.extend_ends(rows, files["fasta"], contig, strand, clen[contig], max_ext=6)
    ident = sum(r["status"] == "identical" for r in rows)
    before, after = M6.flanks(files["fasta"], contig, strand, rows, clen[contig])
    exons = len({r["exon"] for r in rows if r["exon"] != ""})
    return rows, ident, exons, strand, before, after


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--genomedir", required=True)
    ap.add_argument("--curated", required=True, help="step 4 curated_loci.tsv")
    ap.add_argument("--names", help="step 5 names_names.tsv (loci labelled by their new names)")
    ap.add_argument("--species", default="HelSus,HelCha")
    ap.add_argument("--peptides", default=os.path.join(HERE, "8_published_peptides.fa"))
    ap.add_argument("--pad-bp", type=int, default=5000, help="search window around each locus [5000]")
    ap.add_argument("--inside-bp", type=int, default=300,
                    help="the chain's anchor must lie within the locus +/- this [300]")
    ap.add_argument("--min-anchor-bits", type=float, default=15.0)
    ap.add_argument("--min-exon-bits", type=float, default=16.0)
    ap.add_argument("--fill-bits", type=float, default=12.0)
    ap.add_argument("--min-intron-bp", type=int, default=60)
    ap.add_argument("--max-intron-bp", type=int, default=10000)
    ap.add_argument("--frameshift-pen", type=float, default=15.0)
    ap.add_argument("--intron-pen", type=float, default=4.0)
    ap.add_argument("--nosplice-pen", type=float, default=25.0)
    ap.add_argument("--splice-indel-codons", type=int, default=4)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    peps = M6.read_peptides(a.peptides)
    species = [x.strip() for x in a.species.split(",") if x.strip()]
    table = {r["code"]: r for r in N.read_species_table(a.species_table)}
    old2new = {}
    if a.names:
        for r in read_tsv(a.names):
            old2new[r["old_name"]] = r["new_name"]
    loci = defaultdict(list)
    for r in read_tsv(a.curated):
        if r["species"] in species:
            loci[r["species"]].append(dict(r, label=old2new.get(r["name"], r["name"])))
    for sp in species:
        loci[sp].sort(key=lambda r: (r["contig"], int(r["g_start"])))
        N.info(f"{sp}: {len(loci[sp])} curated loci: " + ", ".join(r["label"] for r in loci[sp]))
    tmp = tempfile.mkdtemp(prefix="8a_")
    M_rows, R_rows = [], []
    for sp in species:
        if sp not in table:
            N.warn(f"{sp} not in the species table"); continue
        files = N.find_genome_files(a.genomedir, sp, table[sp]["accession"])
        if not os.path.exists(files["fasta"] + ".fai"):
            M6.run(["samtools", "faidx", files["fasta"]])
        clen = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(files["fasta"] + ".fai")}
        for rec in loci[sp]:
            w0 = max(1, int(rec["g_start"]) - a.pad_bp)
            w1 = min(clen[rec["contig"]], int(rec["g_end"]) + a.pad_bp)
            plus = M6.faidx(files["fasta"], rec["contig"], w0, w1)
            for p in peps:
                out = map_peptide(p["seq"], rec, files, clen, plus, w0, w1, a, tmp)
                base = dict(peptide=p["id"], group=p["group"], label=p["label"], length=len(p["seq"]), species=sp,
                            locus=rec["label"], step4_name=rec["name"], contig=rec["contig"])
                if out is None:
                    M_rows.append(dict(base, identical=0, mapped=0, exons=0, strand="", before="", after="", start="",
                                       end="", note="no hit in the locus"))
                    continue
                rows, ident, ex, strand, before, after = out
                cs = [int(r["codon_start"]) for r in rows if r["codon_start"] != ""]
                ce = [int(r["codon_end"]) for r in rows if r["codon_end"] != ""]
                M_rows.append(dict(base, identical=ident, mapped=len(cs), exons=ex, strand=strand,
                                   before=before or "-", after=after or "-", start=min(cs) if cs else "",
                                   end=max(ce) if ce else "", note=""))
                for r in rows:
                    R_rows.append(dict(peptide=p["id"], species=sp, locus=rec["label"], **r))
    # ranks and best per genome
    best_rows = []
    for p in peps:
        per = {}
        for sp in species:
            ms = sorted((m for m in M_rows if m["peptide"] == p["id"] and m["species"] == sp),
                        key=lambda m: -int(m["identical"]))
            for k, m in enumerate(ms, 1):
                m["rank_in_species"] = k
                m["best_in_species"] = "yes" if k == 1 and int(m["identical"]) > 0 else "no"
            per[sp] = ms[0] if ms else None
        order = sorted((sp for sp in species if per[sp]), key=lambda sp: -int(per[sp]["identical"]))
        top = order[0] if order else ""
        margin = (int(per[order[0]]["identical"]) - int(per[order[1]]["identical"])) if len(order) > 1 else ""
        br = dict(peptide=p["id"], label=p["label"], length=len(p["seq"]), source=p["source"], best_genome=top,
                  margin=margin)
        for sp in species:
            m = per[sp]
            br[f"{sp}_locus"] = m["locus"] if m else "-"
            br[f"{sp}_identical"] = m["identical"] if m else 0
            br[f"{sp}_before"] = m["before"] if m else ""
            br[f"{sp}_after"] = m["after"] if m else ""
        best_rows.append(br)
        N.info(f"{p['label']}: " + "; ".join(
            f"{sp} best {per[sp]['locus']} {per[sp]['identical']}/{len(p['seq'])}" if per[sp] else f"{sp} none"
            for sp in species) + (f" -> {top} (+{margin})" if top else ""))
        for m in M_rows:
            if m["peptide"] == p["id"]:
                m["best_genome"] = top

    def write(fn, cols, rows):
        with open(os.path.join(a.outdir, fn), "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    write("matches.tsv", ["peptide", "group", "label", "length", "species", "locus", "step4_name", "contig", "strand",
                          "start", "end", "identical", "mapped", "exons", "before", "after", "rank_in_species",
                          "best_in_species", "best_genome", "note"], M_rows)
    bcols = ["peptide", "label", "length", "best_genome", "margin"]
    for sp in species:
        bcols += [f"{sp}_locus", f"{sp}_identical", f"{sp}_before", f"{sp}_after"]
    write("best.tsv", bcols + ["source"], best_rows)
    write("residues.tsv", ["peptide", "species", "locus", "qpos", "qres", "exon", "codon_start", "codon_end", "gres",
                           "status", "ins_after", "ext"], R_rows)
    with open(os.path.join(a.outdir, "peptides.tsv"), "w") as fh:
        fh.write("peptide\tgroup\tlabel\tsource\tseq\n")
        for p in peps:
            fh.write(f"{p['id']}\t{p['group']}\t{p['label']}\t{p['source']}\t{p['seq']}\n")
    N.info(f"wrote {a.outdir}/matches.tsv, best.tsv, residues.tsv")


if __name__ == "__main__":
    main()
