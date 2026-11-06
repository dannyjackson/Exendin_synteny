#!/usr/bin/env python3
"""
Build the protein sets for the reciprocal BLAST search (step 4).

FORWARD queries (searched against every genome with tblastn, 4b):
  - every exendin-family protein from the NCBI annotations (2a: ncbi_exendin_aa.fa,
    genes in the interval, near a boundary, or unplaced);
  - every TOGA2 exendin ortholog that passed 3b's filters (tips.tsv, included = yes);
  - the published exendin precursors in 4_reference_peptides.fa (*_precursor records:
    exendin-2, -3, -4). Mature peptides are not used: they are too short to search with.
  Identical sequences are kept once (priority published > NCBI > TOGA2); the other names
  are listed as synonyms. Each query carries a NAME group (lib_exendin_names), which is
  what an ortholog is named after in 4d.
BOUNDARY queries: the NCBI proteins of SLC4A8 and AQP5 from --boundary-species (2a:
  ncbi_boundary_aa.fa), so the boundaries can be found in unannotated genomes too.
REVERSE database (blastx of every forward locus back against it, 4d):
  - one or more whole reference proteomes (--reverse-proteome CODE=protein.faa). These
    supply the competitors: proglucagon, GIP, VIP, PACAP, secretin, GHRH ... A locus
    whose best reverse hit is one of them is not an exendin;
  - all forward queries, plus the published proglucagon precursor.
  Proteome proteins are named by the product in their FASTA header, so a HelSus protein
  called 'exendin-4-like' is in group exendin-4 and is a family member.

Read-only on its inputs. Writes to --outdir:
  forward_queries.faa / query_table.tsv      id species source group family length synonyms description
  boundary_queries.faa
  reverse_db.faa / reverse_table.tsv          same columns as query_table.tsv

Usage
  python3 4a_build_query_set.py \
      --ncbi-dir ../synteny_ncbi/ncbi_loci \
      --toga2-tips ../synteny_toga2/toga2_exendins/tips.tsv \
      --toga2-reference ../synteny_toga2/reference_exendins.tsv \
      --published ${SCRIPTDIR}/4_reference_peptides.fa \
      --reverse-proteome HelSus=${GENOMEDIR}/HelSus_GCF_057233095.1/ncbi_dataset/data/GCF_057233095.1/protein.faa \
      --outdir queries
"""
import argparse, csv, os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

COLS = ["id", "species", "source", "group", "family", "length", "synonyms", "description"]
PRIORITY = {"published": 0, "ncbi": 1, "toga2": 2}


def clean(seq):
    return re.sub(r"[^A-Z]", "", seq.upper().rstrip("*")).replace("*", "X")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncbi-dir", required=True, help="2a output directory")
    ap.add_argument("--toga2-tips", help="tips.tsv from 3b")
    ap.add_argument("--toga2-reference", help="reference_exendins.tsv from 3a (names the TOGA2 tips)")
    ap.add_argument("--toga2-name-map", help="as 3c --name-map")
    ap.add_argument("--published", help="4_reference_peptides.fa")
    ap.add_argument("--reverse-proteome", action="append", default=[], metavar="CODE=FASTA",
                    help="whole proteome for the reverse search; repeat for several")
    ap.add_argument("--boundary-species", default="HelSus",
                    help="species whose boundary proteins are the boundary queries (fallback: any)")
    ap.add_argument("--min-aa", type=int, default=50, help="shorter forward queries are dropped")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    if not a.reverse_proteome:
        N.warn("no --reverse-proteome: the reverse search has no non-exendin competitors, so every forward "
               "locus will look like an exendin. Give at least the HelSus RefSeq proteome.")

    cand = []          # (priority, species, source, group, description, seq, name)
    if a.published:
        for name, (head, seq) in N.read_fasta(a.published).items():
            if "_precursor" not in name:
                continue
            grp = N.name_group(name.replace("EXENDIN", "exendin-").replace("GCG_", "proglucagon "))
            cand.append((0, "published", "published", grp, head.replace(" ", "_"), clean(seq),
                         name.split("|")[0]))
    loci = os.path.join(a.ncbi_dir, "ncbi_loci.tsv")
    pid_row = {}
    with open(loci) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["kind"] == "exendin" and r["protein_id"] != "NA":
                pid_row[(r["species"], r["protein_id"])] = r
    for name, (head, seq) in N.read_fasta(os.path.join(a.ncbi_dir, "ncbi_exendin_aa.fa")).items():
        pid = head.split()[1] if len(head.split()) > 1 else ""
        sp = name.split("_")[0]
        r = pid_row.get((sp, pid))
        grp = r["group"] if r else N.name_group(head)
        cand.append((1, sp, "ncbi", grp, f"{pid} {r['description'] if r else ''} {r['placement'] if r else ''}",
                     clean(seq), name))
    if a.toga2_tips:
        if not a.toga2_reference:
            sys.exit("[error] --toga2-tips needs --toga2-reference (to name the tips)")
        ref_group = N.toga2_reference_groups(a.toga2_reference, a.toga2_name_map)
        with open(a.toga2_tips) as fh:
            for t in csv.DictReader(fh, delimiter="\t"):
                if t["included"] != "yes" or t.get("kind") == "paralog":
                    continue
                cand.append((2, t["query"], "toga2", N.toga2_tip_group(t["ref_genes"], ref_group),
                             f"{t['q_gene']} {t['projection']} {t['status']} {t['orthology_class']}",
                             clean(t["protein"]), t["tip"]))

    # ---- forward queries: family only, deduplicated by sequence ----
    by_seq, dropped = {}, defaultdict(int)
    for c in sorted(cand, key=lambda c: c[0]):
        pr, sp, src, grp, desc, seq, name = c
        if not N.is_family(grp):
            dropped["not a family name: " + grp] += 1
            continue
        if len(seq) < a.min_aa:
            dropped[f"shorter than {a.min_aa} aa"] += 1
            continue
        if seq in by_seq:
            k = by_seq[seq]
            k["synonyms"].append(f"{sp}:{name}")
            if k["group"] != grp:
                N.warn(f"identical sequence named {k['group']} ({k['name']}) and {grp} ({sp}:{name}); "
                       f"kept {k['group']}")
            continue
        by_seq[seq] = dict(name=name, species=sp, source=src, group=grp, description=desc, seq=seq, synonyms=[])
    fw, n = [], defaultdict(int)
    for k in by_seq.values():
        ab = N.group_abbrev(k["group"])
        n[(k["species"], k["source"], ab)] += 1
        k["id"] = f"{k['species']}__{k['source']}__{ab}__{n[(k['species'], k['source'], ab)]}"
        fw.append(k)
    for why, c in dropped.items():
        N.info(f"forward set: dropped {c} protein(s): {why}")

    def write(fasta, table, recs):
        with open(os.path.join(a.outdir, fasta), "w") as fa, open(os.path.join(a.outdir, table), "w") as ft:
            ft.write("\t".join(COLS) + "\n")
            for k in recs:
                fa.write(f">{k['id']}\n{k['seq']}\n")
                ft.write("\t".join([k["id"], k["species"], k["source"], k["group"],
                                    "yes" if N.is_family(k["group"]) else "no", str(len(k["seq"])),
                                    ",".join(k.get("synonyms", [])) or "-",
                                    re.sub(r"\s+", " ", k.get("description", "")).strip() or "-"]) + "\n")
    write("forward_queries.faa", "query_table.tsv", fw)
    by_group = defaultdict(int)
    for k in fw:
        by_group[k["group"]] += 1
    N.info(f"{len(fw)} forward queries: " + ", ".join(f"{g} {c}" for g, c in sorted(by_group.items())))

    # ---- boundary queries ----
    bnd = N.read_fasta(os.path.join(a.ncbi_dir, "ncbi_boundary_aa.fa"))
    chosen = {}
    for name, (head, seq) in bnd.items():
        bname, sp = name[len("boundary_"):].split("__")[0], name.split("__")[1].split("_")[0]
        if bname not in chosen or (sp == a.boundary_species and chosen[bname][0] != a.boundary_species):
            chosen[bname] = (sp, name, seq)
    with open(os.path.join(a.outdir, "boundary_queries.faa"), "w") as fh:
        for bname, (sp, name, seq) in sorted(chosen.items()):
            fh.write(f">{name}\n{clean(seq)}\n")
            N.info(f"boundary query {bname}: {name} ({len(seq)} aa)"
                   + ("" if sp == a.boundary_species else f"  [not from {a.boundary_species}]"))
    if len(chosen) < 2:
        N.warn(f"only {len(chosen)} boundary query/queries: 4d will not place loci between the boundaries")

    # ---- reverse database ----
    rev = []
    for item in a.reverse_proteome:
        code, _, path = item.partition("=")
        ng = defaultdict(int)
        for pid, (head, seq) in N.read_fasta(path).items():
            prod = re.sub(r"\s*\[[^\]]*\]\s*$", "", head[len(pid):]).strip()
            prod = re.sub(r"^(PREDICTED|LOW QUALITY PROTEIN):\s*", "", prod)
            grp = N.name_group(prod)
            ng[grp if not grp.startswith("other") else "other"] += 1
            rev.append(dict(id=f"rp__{code}__{pid}", species=code, source="proteome", group=grp,
                            description=prod, seq=clean(seq)))
        N.info(f"reverse proteome {code}: {sum(ng.values())} proteins; "
               + ", ".join(f"{g} {c}" for g, c in sorted(ng.items()) if g != "other"))
    rev += [dict(k) for k in fw]
    for pr, sp, src, grp, desc, seq, name in cand:
        if src == "published" and not N.is_family(grp):
            rev.append(dict(id=f"published__{N.group_abbrev(grp)}__{name}", species="published",
                            source="published", group=grp, description=desc, seq=seq))
    write("reverse_db.faa", "reverse_table.tsv", rev)
    N.info(f"{len(rev)} proteins in the reverse database")


if __name__ == "__main__":
    main()
