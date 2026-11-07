#!/usr/bin/env python3
"""
Select the curated step-4 genes for a gene tree, and write the files the tree steps need.

Input: the curated loci (4d curated_loci.tsv) and their models (4e models.tsv,
curated_aa.fa, curated_cds.fa). Every curated exendin-family locus in the SLC4A8-AQP5
region is a tip by default: remnants, flagged models and pseudogenes included, since the
question is the whole family in the region. Only these are left out:
  - loci with no 4e model (nothing to align);
  - loci 4d placed 'elsewhere', i.e. outside the region (--placements to change);
  - --exclude NAME,NAME, by hand.
Stop codons are masked so the tree programs accept the sequence: internal '*' -> 'X' in
the protein, stop codons -> NNN in the CDS (the gene stays in; tips.tsv says what was
masked). --exclude-flagged and --exclude-disabled restore the stricter selection, e.g. for
a selection analysis, where pseudogenes and partial models do not belong.
Very short models (single-exon remnants of ~20-60 aa) are kept too; they are mostly gaps
in the alignment and can sit on long, poorly supported branches; read their placement
with that in mind, or leave them out with --exclude.

Outgroup: --outgroup-fasta (protein sequences added as tips, e.g. chicken's NCBI-annotated
GCGL, when the genome model of the outgroup is partial), or else the tips of --outgroup-species (e.g. GalGal: chicken GCGL, the single-copy
ancestral gene of the array outside squamates, as zebra finch GCGL was in the old repo).
With one outgroup sequence the root is placed on that branch by construction: the
outgroup orients the tree, it does not test where the duplications happened. Without one,
root at the midpoint (lib_midpoint_root.py) and treat the root as a display choice.

Usage
  python3 4g_tree_inputs.py --curated reciprocal/curated_loci.tsv --models-dir models \
      --outgroup-species GalGal --outdir tree

Writes to --outdir:
  tree_aa.fa, tree_cds.fa   sequences of the tips (tip = locus name from 4d); --published
                            adds the published exendin-2/-3/-4 precursors (Published_EX2 ...)
                            to tree_aa.fa only
  groups.tsv                tip<TAB>group (name group, as in the synteny plots)
  outgroup.txt              comma-separated outgroup tips (empty: midpoint-root)
  tips.tsv                  every curated locus: included yes/no and why
  name_sources.tsv          (with --reverse-table) tip, code and origin of its name: the
                            protein its name group was taken from, i.e. its best reverse hit
                            in 4d. R = reference RefSeq proteome, N = NCBI annotation of
                            another species, T = TOGA2 projection, P = published precursor,
                            O = outgroup sequence from --outgroup-fasta
"""
import argparse, csv, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--curated", required=True, help="curated_loci.tsv from 4d")
    ap.add_argument("--models-dir", required=True, help="4e output (models.tsv, curated_aa.fa, curated_cds.fa)")
    ap.add_argument("--outgroup-species", default="", help="comma-separated species codes, e.g. GalGal")
    ap.add_argument("--outgroup-fasta",
                    help="protein FASTA of outgroup tips added as they are (e.g. chicken's NCBI GCGL protein), "
                         "used instead of --outgroup-species models; no CDS, so not in tree_cds.fa")
    ap.add_argument("--outgroup-group", default="outgroup", help="colour group of --outgroup-fasta tips")
    ap.add_argument("--published",
                    help="4_reference_peptides.fa: its exendin *_precursor records are added as tips "
                         "(Published_EX2, _EX3, _EX4); protein only, so not in tree_cds.fa")
    ap.add_argument("--reverse-table",
                    help="reverse_table.tsv from 4a: writes name_sources.tsv (where each tip's name came from), "
                         "which 4h --sources prints as a code after each tip")
    ap.add_argument("--placements", default="interval,near_boundary,unplaced",
                    help="4d placements that count as the SLC4A8-AQP5 region (add 'elsewhere' for all)")
    ap.add_argument("--exclude-flagged", action="store_true",
                    help="leave out models whose 4e status is not ok (default: keep every model)")
    ap.add_argument("--exclude-disabled", action="store_true",
                    help="leave out copies with a premature stop (default: keep them, stops masked)")
    ap.add_argument("--exclude", default="", help="comma-separated tip names to leave out")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    with open(a.curated) as fh:
        cur = list(csv.DictReader(fh, delimiter="\t"))
    with open(os.path.join(a.models_dir, "models.tsv")) as fh:
        models = {r["name"]: r for r in csv.DictReader(fh, delimiter="\t")}
    aa = {k: v[1] for k, v in N.read_fasta(os.path.join(a.models_dir, "curated_aa.fa")).items()}
    cdsf = os.path.join(a.models_dir, "curated_cds.fa")
    cds = {k: v[1] for k, v in N.read_fasta(cdsf).items()} if os.path.exists(cdsf) else {}
    excl = set(x for x in a.exclude.split(",") if x)
    og_sp = set(x for x in a.outgroup_species.split(",") if x)

    region = set(x for x in a.placements.split(",") if x)
    rows, keep, masked = [], [], {}
    for r in cur:
        n, why = r["name"], []
        m = models.get(n)
        seq = aa.get(n, "").rstrip("*")
        stopped = "premature_stop" in r.get("state", "") or "*" in seq or (m and "internal_stop" in m["status"])
        if n not in aa or not m:
            why.append("no model (4e)")
        elif not seq:
            why.append("empty model")
        if r.get("placement", "interval") not in region:
            why.append(f"placement {r.get('placement')} (outside the region)")
        if a.exclude_flagged and m and m["status"] != "ok":
            why.append(f"4e status {m['status']} (--exclude-flagged)")
        if a.exclude_disabled and stopped:
            why.append("premature stop (--exclude-disabled)")
        if n in excl:
            why.append("--exclude")
        note = []
        if not why:
            if "*" in seq:
                note.append(f"{seq.count('*')} internal stop(s) masked as X")
            if m and m["status"] != "ok":
                note.append(f"4e status {m['status']}")
            masked[n] = "; ".join(note)
        rows.append([n, r["species"], r["group"], r["class"], m["aa_len"] if m else "0",
                     "yes" if not why else "no", "; ".join(why) or ("; ".join(note) or "-")])
        if not why:
            keep.append(r)

    extra = []
    if a.outgroup_fasta:
        for n, (h, seq) in N.read_fasta(a.outgroup_fasta).items():
            if n in {r["name"] for r in keep}:
                N.warn(f"--outgroup-fasta {n}: same name as a curated tip; skipped"); continue
            extra.append((n, seq.rstrip("*")))
        if og_sp:
            N.info("--outgroup-fasta given: --outgroup-species tips stay in as ordinary tips")
        og_sp = set()
    og = [n for n, _ in extra] + [r["name"] for r in keep if r["species"] in og_sp]
    published = []                    # (tip, seq, group, header)
    if a.published:
        for n, (h, seq) in N.read_fasta(a.published).items():
            base = n.split("|")[0]
            m = re.match(r"EXENDIN(\d+)_precursor$", base)
            if m:
                published.append((f"Published_EX{m.group(1)}", seq.rstrip("*"), f"exendin-{m.group(1)}", h))
        N.info(f"published precursors added as tips: {', '.join(t for t, *_ in published) or 'none found'}")
    if og_sp and not og:
        N.warn(f"no included tip from --outgroup-species {','.join(sorted(og_sp))}: midpoint-root instead")
    with open(os.path.join(a.outdir, "tree_aa.fa"), "w") as fa, \
         open(os.path.join(a.outdir, "tree_cds.fa"), "w") as fc, \
         open(os.path.join(a.outdir, "groups.tsv"), "w") as fg:
        fg.write("tip\tgroup\n")
        for r in keep:
            fa.write(f">{r['name']}\n{aa[r['name']].rstrip('*').replace('*', 'X')}\n")
            if r["name"] in cds:
                c = cds[r["name"]]
                c = c[:len(c) - len(c) % 3]
                cod = [c[i:i + 3] for i in range(0, len(c), 3)]
                if cod and cod[-1].upper() in ("TAA", "TAG", "TGA"):
                    cod = cod[:-1]                      # terminal stop off, as in the protein
                cod = ["NNN" if x.upper() in ("TAA", "TAG", "TGA") else x for x in cod]
                fc.write(f">{r['name']}\n{''.join(cod)}\n")
            fg.write(f"{r['name']}\t{r['group']}\n")
        for n, seq in extra:
            fa.write(f">{n}\n{seq}\n")
            fg.write(f"{n}\t{a.outgroup_group}\n")
        for n, seq, grp, h in published:
            fa.write(f">{n}\n{seq}\n")
            fg.write(f"{n}\t{grp}\n")
    open(os.path.join(a.outdir, "outgroup.txt"), "w").write(",".join(og) + "\n")
    with open(os.path.join(a.outdir, "tips.tsv"), "w") as fh:
        fh.write("name\tspecies\tgroup\tclass\taa_len\tincluded\treason\n")
        for row in rows:
            fh.write("\t".join(map(str, row)) + "\n")

    if a.reverse_table:
        with open(a.reverse_table) as fh:
            rev = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
        code_of = {"proteome": "R", "ncbi": "N", "toga2": "T", "published": "P"}
        with open(os.path.join(a.outdir, "name_sources.tsv"), "w") as fh:
            fh.write("tip\tcode\tsource\tspecies\tprotein\tdescription\n")
            for r in keep:
                x = rev.get(r["rev_best"])
                if x:
                    fh.write("\t".join([r["name"], code_of.get(x["source"], "?"), x["source"], x["species"],
                                        r["rev_best"], x["description"]]) + "\n")
                else:
                    fh.write(f"{r['name']}\t?\tunknown\t-\t{r['rev_best']}\t-\n")
                    N.warn(f"{r['name']}: best reverse hit {r['rev_best']} not in {a.reverse_table}")
            for n, _ in extra:
                fh.write(f"{n}\tO\toutgroup\t-\t{n}\tadded with --outgroup-fasta\n")
            for n, seq, grp, h in published:
                f = h.split("|")
                fh.write(f"{n}\tP\tpublished\t{f[3] if len(f) > 3 else '-'}\t{f[1] if len(f) > 1 else n}\t"
                         f"{f[2] if len(f) > 2 else grp}\n")
        N.info(f"name sources written to {os.path.join(a.outdir, 'name_sources.tsv')}")

    by = {}
    for r in keep:
        by[r["group"]] = by.get(r["group"], 0) + 1
    N.info(f"{len(keep)}/{len(cur)} curated loci in the tree from {len({r['species'] for r in keep})} species: "
           + ", ".join(f"{g} {c}" for g, c in sorted(by.items())))
    for row in rows:
        if row[5] == "no":
            N.info(f"  left out {row[0]}: {row[6]}")
    N.info(f"outgroup: {','.join(og) or 'none (midpoint-root)'}")
    if len(keep) < 4:
        N.warn("fewer than 4 tips: IQ-TREE needs at least 4")


if __name__ == "__main__":
    main()
