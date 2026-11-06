#!/usr/bin/env python3
"""
Build species.tsv: one row per taxon, with the NCBI assembly to download for it.

Taxa come from three places, merged by code (3+3 letters of Genus_species, e.g. HelSus):
  - the tips of the Cactus species tree (--tree), which is the backbone of every plot;
  - the TOGA2 queries (query_* directories of --toga2-run), so step 3 can be drawn for
    every species step 2 draws;
  - --extra "Genus_species,Genus_species" for outgroups outside the tree (drawn below it).

Accessions come from --pinned (1_species_pinned.tsv) when the code is there. Otherwise
`datasets summary genome taxon` is asked for the species, and the assembly picked is:
RefSeq (GCF) with an annotation > any annotated assembly > RefSeq reference/representative
> highest contig N50. Every candidate NCBI returned is written to --candidates so the
choice can be checked. Nothing is downloaded here.

Tree tips that are not Genus_species (e.g. codes such as 'HelSus') need a --pinned row,
which supplies the species name.

Usage
  python3 1a_make_species_table.py \
      --tree ${CACTUSDIR}/Anguimorpha_species_tree_iteration-01_5-23-26.nwk \
      --toga2-run ${TOGADIR}/target_HelSus \
      --pinned ${SCRIPTDIR}/1_species_pinned.tsv \
      --out species.tsv --candidates species_candidates.tsv

Output columns (species.tsv)
  code species accession toga2_query tree_tip source note
  source = pinned | ncbi_auto | none. Edit the file by hand if a choice is wrong;
  every later step reads it.
"""
import argparse, json, os, re, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_exendin_names import TABLE_COLS, info, warn  # noqa: E402


def tips_of(newick):
    s = re.sub(r"\[[^\]]*\]", "", newick)
    return [t.strip().strip("'") for t in re.findall(r"[(,]\s*('[^']+'|[^(),:;]+)", s) if t.strip()]


def code_of(name):
    parts = [p for p in re.split(r"[_\s.]+", name) if p]
    if len(parts) >= 2 and parts[0][:1].isupper():
        return parts[0][:3].capitalize() + parts[1][:3].capitalize()
    return None


def read_pinned(path):
    out = {}
    if not path:
        return out
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if f[0] == "code":
            continue
        out[f[0]] = dict(species=f[1] if len(f) > 1 else "", accession=f[2] if len(f) > 2 else "NA",
                         note=f[3] if len(f) > 3 else "")
    return out


def ncbi_candidates(species):
    if not shutil.which("datasets"):
        return None
    cmd = ["datasets", "summary", "genome", "taxon", species.replace("_", " "), "--as-json-lines"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        warn(f"datasets failed for {species}: {r.stderr.strip()[:200]}")
        return []
    out = []
    for line in r.stdout.splitlines():
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            continue
        ai, st = d.get("assembly_info", {}), d.get("assembly_stats", {})
        org = d.get("organism", {}).get("organism_name", "")
        out.append(dict(accession=d.get("accession", ""), organism=org,
                        annotated="annotation_info" in d,
                        refseq=d.get("accession", "").startswith("GCF_"),
                        category=ai.get("refseq_category", ""), level=ai.get("assembly_level", ""),
                        name=ai.get("assembly_name", ""), contig_n50=int(st.get("contig_n50", 0) or 0),
                        released=ai.get("release_date", "")))
    # an exact species match only (datasets also returns subspecies / other taxa)
    want = species.replace("_", " ").lower()
    exact = [c for c in out if c["organism"].lower().startswith(want)]
    return exact or out


def rank(c):
    return (not (c["refseq"] and c["annotated"]), not c["annotated"],
            not (c["refseq"] and c["category"] in ("reference genome", "representative genome")),
            -c["contig_n50"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True, help="Cactus species tree (Newick)")
    ap.add_argument("--toga2-run", help="${TOGADIR}/target_HelSus (its query_* directories)")
    ap.add_argument("--pinned", help="1_species_pinned.tsv")
    ap.add_argument("--extra", default="", help="comma-separated Genus_species outside the tree")
    ap.add_argument("--out", required=True)
    ap.add_argument("--candidates", help="TSV of every NCBI candidate, for checking the choice")
    a = ap.parse_args()

    pinned = read_pinned(a.pinned)
    rows = {}

    def add(code, species, **kw):
        r = rows.setdefault(code, {c: "" for c in TABLE_COLS})
        r["code"] = code
        r["species"] = r["species"] or species
        for k, v in kw.items():
            if v and not r[k]:
                r[k] = v

    for tip in tips_of(open(a.tree).read()):
        c = code_of(tip)
        if c is None:
            if tip in pinned:
                add(tip, pinned[tip]["species"], tree_tip=tip)
            else:
                warn(f"tree tip '{tip}' is not Genus_species and has no --pinned row: skipped")
            continue
        add(c, tip if "_" in tip else tip.replace(" ", "_"), tree_tip=tip)
    if a.toga2_run:
        for d in sorted(os.listdir(a.toga2_run)):
            if d.startswith("query_"):
                q = d[6:]
                if q in rows:
                    rows[q]["toga2_query"] = q
                else:
                    add(q, pinned.get(q, {}).get("species", ""), toga2_query=q,
                        note="TOGA2 query not in the tree")
    for sp in [s for s in a.extra.split(",") if s.strip()]:
        c = code_of(sp.strip())
        if c:
            add(c, sp.strip(), note="outgroup, not in the tree")

    cand_rows = []
    for code, r in sorted(rows.items()):
        if code in pinned:
            p = pinned[code]
            r["accession"], r["source"] = p["accession"], "pinned"
            r["species"] = r["species"] or p["species"]
            r["note"] = "; ".join(x for x in (r["note"], p["note"]) if x)
            continue
        if not r["species"]:
            r["accession"], r["source"] = "NA", "none"
            r["note"] = "; ".join(x for x in (r["note"], "species name unknown: add a --pinned row") if x)
            continue
        cands = ncbi_candidates(r["species"])
        if cands is None:
            r["accession"], r["source"] = "NA", "none"
            r["note"] = "; ".join(x for x in (r["note"], "datasets not on PATH") if x)
            continue
        if not cands:
            r["accession"], r["source"] = "NA", "none"
            r["note"] = "; ".join(x for x in (r["note"], "no NCBI assembly") if x)
            continue
        cands.sort(key=rank)
        best = cands[0]
        r["accession"], r["source"] = best["accession"], "ncbi_auto"
        r["note"] = "; ".join(x for x in (r["note"], f"{best['name']} {best['level']}, "
                              f"{'annotated' if best['annotated'] else 'no annotation'}, "
                              f"contig N50 {best['contig_n50']:,}; {len(cands)} candidates") if x)
        for c in cands:
            cand_rows.append([code, r["species"], c["accession"], c["name"], c["organism"],
                              "yes" if c["annotated"] else "no", c["category"], c["level"],
                              c["contig_n50"], c["released"], "chosen" if c is best else ""])
        if r["toga2_query"]:
            warn(f"{code} is a TOGA2 query but not pinned: check that {best['accession']} is the assembly "
                 f"TOGA2 used, or contig names will not match between steps 2 and 3")

    with open(a.out, "w") as fh:
        fh.write("\t".join(TABLE_COLS) + "\n")
        for code in sorted(rows):
            fh.write("\t".join(rows[code][c] for c in TABLE_COLS) + "\n")
    if a.candidates:
        with open(a.candidates, "w") as fh:
            fh.write("code\tspecies\taccession\tassembly_name\torganism\tannotated\trefseq_category\t"
                     "level\tcontig_n50\treleased\tchosen\n")
            for row in cand_rows:
                fh.write("\t".join(map(str, row)) + "\n")
    n_na = sum(r["accession"] in ("", "NA") for r in rows.values())
    info(f"{len(rows)} taxa written to {a.out}; {n_na} without an accession")
    for code in sorted(rows):
        r = rows[code]
        info(f"  {code:8s} {r['species']:32s} {r['accession']:18s} {r['source']:9s} "
             f"{'TOGA2 ' if r['toga2_query'] else '      '}{'tree' if r['tree_tip'] else '    '}")


if __name__ == "__main__":
    main()
