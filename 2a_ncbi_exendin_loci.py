#!/usr/bin/env python3
"""
Find, in each species' NCBI annotation, the exendin-family genes and the two boundary
genes of the conserved interval (SLC4A8 ... AQP5), and write one table for the step-2
synteny plot (2b_plot_ncbi_synteny.py) and the step-4 query set.

Inputs
  - species.tsv (1a) and ${GENOMEDIR} (1b): genomic.gff and protein.faa of each species.
    A species with no GFF (assembly not annotated by NCBI) gets one 'no_annotation' row.
  - --pattern: a gene is an exendin-family candidate if its symbol, description or any
    transcript/CDS product matches (case-insensitive), and --exclude does not.
    Default 'exendin|glucagon' minus 'receptor', so genes the annotation calls
    glucagon-like (GCGL, the ancestral gene of the array) are kept. Real proglucagon
    also matches; it is coloured 'proglucagon' and normally sits elsewhere in the genome.

Boundary genes. Each is a gene whose symbol equals NAME (exact) or whose
symbol/description matches REGEX (e.g. 'solute carrier family 4 member 8-like').
Among candidates, the SLC4A8-AQP5 pair on one contig with the most candidate
exendin-family genes between them is used (then exact names, then the shortest
interval). In a fragmented assembly with no such pair, each boundary is chosen
separately, by the number of family genes within --max-dist-kb.

Placement of every family gene
  interval        on the boundary-pair contig, between the boundaries (+/- --pad-kb)
  near_boundary   on a contig with only one boundary, within --max-dist-kb of it
  unplaced        exendin/GCGL-named gene on a contig with no boundary gene
  elsewhere       anything else (proglucagon on its own chromosome, distant copies)
2b draws the first three by default.

Colour group: lib_exendin_names.name_group() of the description/product, so
'exendin-2' and 'exendin-2-like' are one group. Label: 'Ex2-like' etc.

Read-only on ${GENOMEDIR}. Writes to --outdir:
  ncbi_loci.tsv          one row per drawn gene (family and boundary) and per species
                         with no annotation
  ncbi_exendin_aa.fa     longest protein of every family gene (interval, near_boundary,
                         unplaced), named <code>_<group>_<symbol>
  ncbi_boundary_aa.fa    longest protein of each chosen boundary gene, named
                         boundary_<NAME>__<code>_<symbol>

Usage
  python3 2a_ncbi_exendin_loci.py --species-table species.tsv --genomedir ${GENOMEDIR} \
      --outdir ncbi_loci 2> ncbi_loci.log
"""
import argparse, gzip, os, re, sys
from collections import defaultdict
from urllib.parse import unquote

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_exendin_names import (name_group, short_label, group_abbrev, is_family, read_species_table,  # noqa: E402
                               find_genome_files, read_fasta, info, warn)

DEFAULT_BOUNDARIES = "SLC4A8=SLC4A8|solute carrier family 4 member 8,AQP5=AQP5|aquaporin[- ]?5"
COLS = ["species", "kind", "gene_id", "symbol", "description", "group", "label", "contig", "start", "end",
        "strand", "biotype", "protein_id", "aa_len", "placement", "note"]


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def parse_gff(path):
    """gene/pseudogene records with symbol, description, products and protein IDs of their CDS."""
    genes, parent = {}, {}
    prods, prots = defaultdict(set), defaultdict(set)
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 9:
                continue
            at = {k: unquote(v) for k, v in (kv.split("=", 1) for kv in t[8].strip().strip(";").split(";")
                                              if "=" in kv)}
            fid, par = at.get("ID", ""), at.get("Parent", "").split(",")[0]
            if t[2] in ("gene", "pseudogene"):
                genes[fid] = dict(gene_id=fid, contig=t[0], start=int(t[3]), end=int(t[4]), strand=t[6],
                                  symbol=at.get("gene") or at.get("Name") or fid,
                                  description=at.get("description", ""),
                                  biotype=at.get("gene_biotype") or ("pseudogene" if t[2] == "pseudogene"
                                                                     else "NA"))
            else:
                if fid and par:
                    parent[fid] = par
                key = par or fid
                if "product" in at:
                    prods[key].add(at["product"])
                if t[2] == "CDS" and "protein_id" in at:
                    prots[key].add(at["protein_id"])

    def top(x):
        seen = 0
        while x not in genes and x in parent and seen < 5:
            x, seen = parent[x], seen + 1
        return x
    for d, src in ((prods, "products"), (prots, "proteins")):
        for k, v in d.items():
            g = top(k)
            if g in genes:
                genes[g].setdefault(src, set()).update(v)
    for g in genes.values():
        g.setdefault("products", set()); g.setdefault("proteins", set())
    return genes


def text_of(g):
    return " | ".join([g["symbol"], g["description"]] + sorted(g["products"]))


def best_protein(g, prot):
    best = None
    for pid in g["proteins"]:
        if pid in prot and (best is None or len(prot[pid][1]) > len(prot[best][1])):
            best = pid
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-table", required=True, help="species.tsv from 1a")
    ap.add_argument("--genomedir", required=True)
    ap.add_argument("--pattern", default="exendin|glucagon", help="regex for family genes")
    ap.add_argument("--exclude", default="receptor", help="regex; matching genes are not family genes")
    ap.add_argument("--boundaries", default=DEFAULT_BOUNDARIES, help="two NAME=REGEX, left one first")
    ap.add_argument("--pad-kb", type=float, default=500.0,
                    help="family genes this far outside the boundary pair still count as interval "
                         "(HelSus has copies 118 and 379 kb from ANC)")
    ap.add_argument("--max-dist-kb", type=float, default=1000.0,
                    help="for a contig with one boundary: how far from it a family gene is near_boundary")
    ap.add_argument("--max-interval-kb", type=float, default=5000.0,
                    help="boundary pairs further apart than this are not considered")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    pat, excl = re.compile(a.pattern, re.I), re.compile(a.exclude, re.I) if a.exclude else None
    bspec = []
    for item in a.boundaries.split(","):
        name, _, rx = item.partition("=")
        bspec.append((name, re.compile(rf"(?<![A-Za-z0-9])(?:{rx or re.escape(name)})(?![A-Za-z0-9])", re.I)))
    if len(bspec) != 2:
        sys.exit("[error] --boundaries needs exactly two genes")
    pad, maxd, maxi = a.pad_kb * 1000, a.max_dist_kb * 1000, a.max_interval_kb * 1000
    os.makedirs(a.outdir, exist_ok=True)

    rows, fam_fa, bnd_fa = [], [], []
    for sp in read_species_table(a.species_table):
        code = sp["code"]
        if sp["accession"] in ("", "NA"):
            rows.append(dict(species=code, kind="none", placement="no_genome", note="no accession in species.tsv"))
            continue
        files = find_genome_files(a.genomedir, code, sp["accession"])
        if not files["gff"]:
            rows.append(dict(species=code, kind="none", placement="no_annotation",
                             note=f"{sp['accession']}: no NCBI GFF"))
            info(f"{code}: no NCBI annotation ({sp['accession']})")
            continue
        genes = parse_gff(files["gff"])
        prot = read_fasta(files["protein"]) if files["protein"] else {}
        if not files["protein"]:
            warn(f"{code}: no protein.faa; proteins will be missing from the step-4 query set")

        fam = [g for g in genes.values()
               if pat.search(text_of(g)) and not (excl and excl.search(text_of(g)))]
        for g in fam:
            txt = g["description"] or next(iter(sorted(g["products"])), "") or g["symbol"]
            g["group"], g["label"] = name_group(txt), short_label(txt)
        famfam = [g for g in fam if is_family(g["group"])]

        cands = {}
        for name, rx in bspec:
            cands[name] = []
            for g in genes.values():
                exact = g["symbol"].upper() == name.upper()
                if exact or rx.search(g["symbol"]) or rx.search(g["description"]):
                    if g["biotype"] in ("protein_coding", "NA") or exact:
                        cands[name].append(dict(g, exact=exact))
        (n0, _), (n1, _) = bspec

        def n_between(c, lo, hi):
            return sum(1 for g in famfam if g["contig"] == c and g["end"] >= lo and g["start"] <= hi)

        def near_count(b):
            return sum(1 for g in famfam if g["contig"] == b["contig"]
                       and min(abs(g["start"] - b["end"]), abs(b["start"] - g["end"])) <= maxd)

        pairs = []
        for s in cands[n0]:
            for q in cands[n1]:
                if s["contig"] != q["contig"]:
                    continue
                lo, hi = min(s["start"], q["start"]), max(s["end"], q["end"])
                if hi - lo > maxi:
                    continue
                pairs.append(((-n_between(s["contig"], lo, hi), not (s["exact"] and q["exact"]), hi - lo), s, q))
        chosen = {}
        if pairs:
            _, s, q = min(pairs, key=lambda p: p[0])
            chosen = {n0: s, n1: q}
            info(f"{code}: boundary pair on {s['contig']}: {n0}={s['symbol']} ({s['start']}-{s['end']}), "
                 f"{n1}={q['symbol']} ({q['start']}-{q['end']}); {-min(pairs, key=lambda p: p[0])[0][0]} "
                 f"family gene(s) between; {len(pairs)} candidate pair(s)")
        else:
            for name, _ in bspec:
                if cands[name]:
                    chosen[name] = min(cands[name], key=lambda b: (-near_count(b), not b["exact"]))
                    b = chosen[name]
                    info(f"{code}: {name} = {b['symbol']} {b['contig']}:{b['start']}-{b['end']} "
                         f"(no pair on one contig; {near_count(b)} family gene(s) within {a.max_dist_kb:g} kb)")
                else:
                    info(f"{code}: no {name} candidate in the annotation")
        for name, b in chosen.items():
            if not b["exact"]:
                warn(f"{code}: {name} chosen by description, not symbol: {b['symbol']} '{b['description']}'")

        bcontigs = defaultdict(list)
        for name, b in chosen.items():
            bcontigs[b["contig"]].append(b)
        pair_ctg = chosen[n0]["contig"] if pairs else None
        if pairs:
            lo = min(b["start"] for b in chosen.values()) - pad
            hi = max(b["end"] for b in chosen.values()) + pad

        counts = defaultdict(int)
        for g in sorted(fam, key=lambda g: (g["contig"], g["start"])):
            if g["contig"] == pair_ctg and g["end"] >= lo and g["start"] <= hi:
                pl = "interval"
            elif g["contig"] in bcontigs and g["contig"] != pair_ctg and any(
                    min(abs(g["start"] - b["end"]), abs(b["start"] - g["end"])) <= maxd for b in bcontigs[g["contig"]]):
                pl = "near_boundary"
            elif g["contig"] not in bcontigs and is_family(g["group"]):
                pl = "unplaced"
            else:
                pl = "elsewhere"
            counts[pl] += 1
            pid = best_protein(g, prot)
            rows.append(dict(species=code, kind="exendin", gene_id=g["gene_id"], symbol=g["symbol"],
                             description=g["description"] or "; ".join(sorted(g["products"]))[:80],
                             group=g["group"], label=g["label"], contig=g["contig"], start=g["start"],
                             end=g["end"], strand=g["strand"], biotype=g["biotype"], protein_id=pid or "NA",
                             aa_len=len(prot[pid][1]) if pid else 0, placement=pl, note=""))
            if pid and pl != "elsewhere" and is_family(g["group"]):
                fam_fa.append((f"{code}_{group_abbrev(g['group'])}_{g['symbol']}",
                               f"{pid} {g['description']} [{code} {g['contig']}:{g['start']}-{g['end']}]",
                               prot[pid][1]))
        for name, b in chosen.items():
            pid = best_protein(b, prot)
            rows.append(dict(species=code, kind="boundary", gene_id=b["gene_id"], symbol=b["symbol"],
                             description=b["description"], group=name, label=name, contig=b["contig"],
                             start=b["start"], end=b["end"], strand=b["strand"], biotype=b["biotype"],
                             protein_id=pid or "NA", aa_len=len(prot[pid][1]) if pid else 0,
                             placement="boundary_pair" if pairs else "boundary_single",
                             note="exact symbol" if b["exact"] else "matched by description"))
            if pid:
                bnd_fa.append((f"boundary_{name}__{code}_{b['symbol']}", f"{pid} {b['description']}", prot[pid][1]))
        info(f"{code}: {len(fam)} gene(s) match /{a.pattern}/: "
             + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
        for g in fam:
            if g["group"].startswith("other"):
                warn(f"{code}: {g['symbol']} matched /{a.pattern}/ but its name gives no family group: "
                     f"'{g['description']}'")

    with open(os.path.join(a.outdir, "ncbi_loci.tsv"), "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in COLS) + "\n")
    for fn, recs in (("ncbi_exendin_aa.fa", fam_fa), ("ncbi_boundary_aa.fa", bnd_fa)):
        with open(os.path.join(a.outdir, fn), "w") as fh:
            used = defaultdict(int)
            for n, h, s in recs:
                used[n] += 1
                nn = n if used[n] == 1 else f"{n}_{used[n]}"
                fh.write(f">{nn} {h}\n{s.rstrip('*')}\n")
    info(f"wrote {len(rows)} rows, {len(fam_fa)} family proteins, {len(bnd_fa)} boundary proteins to {a.outdir}")


if __name__ == "__main__":
    main()
