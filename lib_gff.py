#!/usr/bin/env python3
"""
Minimal GFF3 reader for step 6. Not run on its own.

  load_region(gff, contig, lo, hi)  genes overlapping contig:lo-hi, each with its
                                    transcripts, exons and CDS (NCBI RefSeq GFF3 layout:
                                    gene/pseudogene -> mRNA/transcript -> exon/CDS)
  gene_label(g)                     'LOC146345366 (uncharacterized protein ...)' style label
  overlaps(genes, lo, hi)           genes (and the exon/CDS numbers) that overlap lo-hi
  find_genes(gff, specs)            genes named by specs ('SLC4A8', 'AQP5=LOC_00012326'),
                                    anywhere in the file: {name: (contig, start, end, strand, id)}
  sequence_regions(gff)             {contig: length} from ##sequence-region lines

Works on NCBI RefSeq GFF3 and on the EviAnn GFF3 TOGA2 projects from (gene names there
are in Note= or Name=; IDs like LOC_00012325).
"""
import gzip, os
from urllib.parse import unquote

GENE_TYPES = {"gene", "pseudogene"}
TX_TYPES = {"mRNA", "transcript", "primary_transcript", "lnc_RNA", "ncRNA", "tRNA", "rRNA",
            "V_gene_segment", "C_gene_segment"}


def _attrs(s):
    out = {}
    for kv in s.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = unquote(v)
    return out


def load_region(gff, contig, lo, hi):
    """Genes overlapping contig:lo-hi -> list of dicts
    {id, name, description, biotype, start, end, strand, tx: [{id, product, exons, cds}]}.
    Exons and CDS are lists of (start, end), sorted. Transcripts without exons get their
    CDS (or their own span) as exons."""
    op = gzip.open(gff, "rt") if gff.endswith(".gz") else open(gff)
    genes, tx, parts = {}, {}, []
    with op as fh:
        for line in fh:
            if line.startswith("#") or not line.startswith(contig + "\t"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != contig:
                continue
            s, e = int(f[3]), int(f[4])
            if e < lo or s > hi:
                continue
            a = _attrs(f[8])
            if f[2] in GENE_TYPES:
                genes[a.get("ID", f"g{len(genes)}")] = dict(
                    id=a.get("ID", ""), name=a.get("Name") or a.get("gene") or a.get("ID", "?"),
                    description=a.get("description") or a.get("Note") or a.get("product", ""),
                    biotype=a.get("gene_biotype", f[2]),
                    start=s, end=e, strand=f[6], tx=[])
            elif f[2] in TX_TYPES:
                tx[a.get("ID", f"t{len(tx)}")] = dict(id=a.get("ID", ""), parent=a.get("Parent", ""),
                                                     product=a.get("product") or a.get("Note", ""), start=s, end=e,
                                                     exons=[], cds=[])
            elif f[2] in ("exon", "CDS"):
                parts.append((f[2], a.get("Parent", ""), s, e, a.get("product", "")))
    for kind, par, s, e, prod in parts:
        for p in par.split(","):
            if p in tx:
                tx[p]["exons" if kind == "exon" else "cds"].append((s, e))
                if prod and not tx[p]["product"]:
                    tx[p]["product"] = prod
            elif p in genes and kind == "exon":         # exon directly under a gene (some ncRNA/pseudogenes)
                t = tx.setdefault(p + "_tx", dict(id=p + "_tx", parent=p, product="", start=s, end=e,
                                                  exons=[], cds=[]))
                t["exons"].append((s, e))
    for t in tx.values():
        t["exons"] = sorted(set(t["exons"])) or sorted(set(t["cds"])) or [(t["start"], t["end"])]
        t["cds"] = sorted(set(t["cds"]))
        g = genes.get(t["parent"])
        if g is not None:
            g["tx"].append(t)
    out = []
    for g in genes.values():
        if not g["tx"]:
            g["tx"] = [dict(id=g["id"], parent=g["id"], product="", start=g["start"], end=g["end"],
                            exons=[(g["start"], g["end"])], cds=[])]
        # representative transcript: the one with the longest CDS, else the longest
        g["tx"].sort(key=lambda t: (-sum(e - s + 1 for s, e in t["cds"]),
                                    -sum(e - s + 1 for s, e in t["exons"])))
        if not g["description"]:
            g["description"] = g["tx"][0]["product"]
        out.append(g)
    return sorted(out, key=lambda g: g["start"])


def gene_label(g, width=48):
    d = g["description"] or g["tx"][0]["product"]
    if d and d.lower() != g["name"].lower():
        d = d if len(d) <= width else d[:width - 1] + "…"
        return f"{g['name']} ({d})"
    return g["name"]


def overlaps(genes, lo, hi):
    """[(gene, 'exon 2 of 3; CDS 2 of 3' | 'intron' | ...)] for genes overlapping lo-hi."""
    out = []
    for g in genes:
        if g["end"] < lo or g["start"] > hi:
            continue
        t = g["tx"][0]
        n = len(t["exons"])
        order = list(range(n)) if g["strand"] != "-" else list(range(n - 1, -1, -1))
        ex = [k + 1 for k, i in enumerate(order) if t["exons"][i][1] >= lo and t["exons"][i][0] <= hi]
        cd = any(e >= lo and s <= hi for s, e in t["cds"])
        if ex:
            where = f"exon {','.join(map(str, ex))} of {n}" + (" (coding)" if cd else " (non-coding)")
        else:
            where = "intron"
        out.append((g, where))
    return out


def _names(a):
    """Every name a GFF record goes by (upper case)."""
    out = set()
    for k in ("ID", "Name", "gene", "gene_name", "Note", "description"):
        v = a.get(k, "")
        if not v:
            continue
        out.add(v.upper())
        if k in ("ID", "Name"):
            out.add(v.split(":")[-1].upper())        # 'gene:LOC_00012326'
        if k in ("Note", "description"):             # 'similar to SLC4A8: ...'
            out.update(w.strip(",;:()").upper() for w in v.replace("=", " ").split())
    return out


def find_genes(gff, specs):
    """specs: ['SLC4A8', 'AQP5=LOC_00012326'] -> {spec_name: (contig, start, end, strand, id)}.
    'NAME=ID' looks the gene up by ID; a bare NAME by Name/gene/ID, then by a word of Note.
    A name matching several genes keeps the first exact Name/gene match."""
    want = []
    for sp in specs:
        n, _, gid = sp.partition("=")
        want.append((n.strip(), gid.strip()))
    found, weak = {}, {}
    op = gzip.open(gff, "rt") if gff.endswith(".gz") else open(gff)
    with op as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in GENE_TYPES:
                continue
            a = _attrs(f[8])
            rec = (f[0], int(f[3]), int(f[4]), f[6], a.get("ID", ""))
            ids = {a.get("ID", "").upper(), a.get("ID", "").split(":")[-1].upper()}
            strong = {a.get(k, "").upper() for k in ("Name", "gene", "gene_name")} | ids
            names = _names(a)
            for n, gid in want:
                if gid:
                    if gid.upper() in ids:
                        found[n] = rec
                elif n.upper() in strong and n not in found:
                    found[n] = rec
                elif n.upper() in names and n not in weak:
                    weak[n] = rec
    for n, _ in want:
        if n not in found and n in weak:
            found[n] = weak[n]
    return found


def sequence_regions(gff):
    out = {}
    op = gzip.open(gff, "rt") if gff.endswith(".gz") else open(gff)
    with op as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                f = line.split()
                if len(f) >= 4:
                    out[f[1]] = int(f[3])
            elif not line.startswith("#"):
                break
    return out


def ref_gene_notes(gff):
    """{gene ID: description} of a reference GFF (EviAnn: Note=), IDs as TOGA2 uses them."""
    out = {}
    op = gzip.open(gff, "rt") if gff.endswith(".gz") else open(gff)
    with op as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in GENE_TYPES:
                continue
            a = _attrs(f[8])
            gid = a.get("ID", "")
            d = a.get("Name") if a.get("Name") and a.get("Name") != gid else ""
            note = a.get("description") or a.get("Note") or ""
            out[gid] = out[gid.split(":")[-1]] = (f"{d}: {note}" if d and note else d or note)
    return out


def load_toga2_query(qdir, notes=None, merge_bp=20000):
    """TOGA2 projections of one query (query_<code>/query_annotation.bed, BED12) as genes
    {contig: [gene]}, in load_region's format. A projection belongs to the reference gene
    in its name (LOC_00012325-mRNA-1.xxx -> LOC_00012325); projections of one reference gene
    on one contig within merge_bp of each other are one gene, named after the reference
    gene, with its reference description (notes) and TOGA2 status (loss_summary.tsv)."""
    def op(p):
        return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

    def first(name):
        for p in (os.path.join(qdir, name), os.path.join(qdir, name + ".gz")):
            if os.path.exists(p):
                return p
        return None
    status = {}
    p = first("loss_summary.tsv")
    if p:
        with op(p) as fh:
            for line in fh:
                t = line.rstrip("\n").split("\t")
                if len(t) >= 3 and t[0] == "PROJECTION":
                    status[t[1]] = t[2]
    p = first("query_annotation.bed")
    if not p:
        raise FileNotFoundError(f"no query_annotation.bed in {qdir}")
    rank = {"FI": 0, "I": 1, "PI": 2, "UL": 3, "M": 4, "L": 5, "PG": 6, "N": 7}
    by = {}
    with op(p) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) < 12:
                continue
            c, s0, e0, name, strand = t[0], int(t[1]), int(t[2]), t[3], t[5]
            ts, te = int(t[6]), int(t[7])
            sizes = [int(x) for x in t[10].strip(",").split(",") if x]
            starts = [int(x) for x in t[11].strip(",").split(",") if x]
            exons = [(s0 + b + 1, s0 + b + z) for b, z in zip(starts, sizes)]
            cds = [(max(s, ts + 1), min(e, te)) for s, e in exons if min(e, te) >= max(s, ts + 1)]
            ref = name.split("-mRNA")[0].split("#")[0]
            st = status.get(name, "NA")
            tx = dict(id=name, parent=ref, product="", start=s0 + 1, end=e0, exons=exons, cds=cds, status=st)
            key = (c, ref, strand)
            for g in by.setdefault(key, []):
                if s0 + 1 <= g["end"] + merge_bp and e0 >= g["start"] - merge_bp:
                    g["tx"].append(tx); g["start"] = min(g["start"], s0 + 1); g["end"] = max(g["end"], e0)
                    break
            else:
                by[key].append(dict(id=f"{ref}@{c}:{s0 + 1}", name=ref, description="", biotype="toga2",
                                    start=s0 + 1, end=e0, strand=strand, tx=[tx], contig=c))
    out = {}
    for (c, ref, _), gs in by.items():
        for g in gs:
            g["tx"].sort(key=lambda x: (rank.get(x["status"], 9), -sum(e - s + 1 for s, e in x["cds"])))
            best = g["tx"][0]["status"]
            d = (notes or {}).get(ref, "")
            g["description"] = (d + " " if d else "") + f"[TOGA2 {best}]"
            g["status"] = best
            out.setdefault(c, []).append(g)
    for c in out:
        out[c].sort(key=lambda g: g["start"])
    return out
