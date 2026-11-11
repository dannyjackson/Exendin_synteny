#!/usr/bin/env python3
"""
One table of every exendin-family gene in one genome (default HelSus): its new name
(10f), position, TOGA2 facts, sequence, and which published exendin its mature peptide
is (10_toga2_exendin_tree.md, step 6).

Joined from
  names/names.tsv            10f: new name, old tip, class, naming rule, coordinates
  collect/candidates.tsv     10a: TOGA2 run, other runs that found it, ortholog/paralog,
                             status, reference gene, projection, strand
  curated/exendin_log.tsv    10b: curation actions (trim, added model ...)
  curated/exendin_aa.fa      10b: the protein that went into the tree
and the mature peptide, its best published match and a verdict, from the code of
5a_assign_nomenclature.py against --reference (4_reference_peptides.fa).

Rows: GCGL, Exendin_anc and the lettered copies of --species (GCG and GIP are left out),
in chromosome order. Writes --out (TSV) and prints a short version.

Usage
  python3 10g_gene_table.py --species HelSus --names names/names.tsv \
      --candidates collect/candidates.tsv --curation-log curated/exendin_log.tsv \
      --aa curated/exendin_aa.fa --reference ${SCRIPTDIR}/4_reference_peptides.fa \
      --out names/HelSus_exendins.tsv
"""
import argparse, csv, importlib.util, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("nomen", os.path.join(HERE, "5a_assign_nomenclature.py"))
NM = importlib.util.module_from_spec(spec)
spec.loader.exec_module(NM)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--species", default="HelSus")
    ap.add_argument("--names", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--curation-log")
    ap.add_argument("--aa", required=True, help="curated protein FASTA (tree tip names)")
    ap.add_argument("--reference", required=True, help="4_reference_peptides.fa (published *_mature peptides)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    names = [r for r in read_tsv(a.names) if r["species"] == a.species and r["class"] not in ("GCG", "GIP")]
    if not names:
        sys.exit(f"[error] no {a.species} genes in {a.names}")
    cand = {r["tip"]: r for r in read_tsv(a.candidates)}
    log = {}
    if a.curation_log:
        for r in read_tsv(a.curation_log):
            log[r["new_name"]] = r
    aa = {h.split()[0]: s.replace("*", "X") for h, s in NM.read_fasta(a.aa)}
    refs = [(h.split("|")[0], s) for h, s in NM.read_fasta(a.reference) if h.split("|")[0].endswith("_mature")]

    rows = []
    for r in names:
        t = r["old_tip"]
        lg = log.get(t, {})
        c = cand.get(lg.get("tip", t), {})
        prot = aa.get(t, "")
        pep, how, best, ident, verdict = "", "", "", "", ""
        if prot and refs:
            pep, how = NM.mature(prot, refs)
            scored = sorted(((max(NM.compare(pep, s), NM.compare(pep[:len(s)], s), key=lambda x: x[0])[0], n)
                             for n, s in refs), reverse=True)
            ident, best = scored[0]
            base = best.replace("_mature", "")
            verdict = (f"published {base} exactly" if ident >= 100 else
                       f"{base}-type ({ident:.0f}%)" if ident >= 80 else f"no published match (closest {base}, {ident:.0f}%)")
        start, end = r.get("start") or c.get("start", ""), r.get("end") or c.get("end", "")
        rows.append(dict(
            name=r["new_name"], old_tip=t, cls=r["class"], position_rule=r["rule"],
            contig=r.get("contig") or c.get("contig", ""), start=start, end=end, strand=c.get("strand", ""),
            locus_kb=f"{(int(end) - int(start) + 1) / 1000:.1f}" if start and end else "",
            source=lg.get("origin", "TOGA2") if lg.get("origin", "TOGA2") != "TOGA2" else f"TOGA2 {c.get('run', '')}".strip(),
            also_in=c.get("also_in", ""), kind=c.get("kind", ""), status=c.get("status", ""),
            ref_gene=c.get("ref_gene", ""), projection=c.get("projection", ""),
            curation=lg.get("actions", ""), aa_len=len(prot), internal_stops=prot.count("X"),
            mature_peptide=pep, peptide_found_by=how, best_published=best.replace("_mature", ""),
            identity=f"{ident:.1f}" if ident != "" else "", verdict=verdict, protein=prot))
    rows.sort(key=lambda x: (x["contig"], int(x["start"]) if str(x["start"]).isdigit() else 0))
    cols = list(rows[0])
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for x in rows:
            fh.write("\t".join(str(x[c]) for c in cols) + "\n")
    show = ["name", "contig", "start", "end", "strand", "status", "aa_len", "verdict", "mature_peptide"]
    w = {c: max(len(c), *(len(str(x[c])) for x in rows)) for c in show}
    print("  ".join(c.ljust(w[c]) for c in show))
    for x in rows:
        print("  ".join(str(x[c]).ljust(w[c]) for c in show))
    print(f"[info] {len(rows)} {a.species} genes; full table (TOGA2 projection, reference gene, curation, protein) "
          f"in {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
