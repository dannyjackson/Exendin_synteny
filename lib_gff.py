#!/usr/bin/env python3
"""
Minimal GFF3 reader for step 6. Not run on its own.

  load_region(gff, contig, lo, hi)  genes overlapping contig:lo-hi, each with its
                                    transcripts, exons and CDS (NCBI RefSeq GFF3 layout:
                                    gene/pseudogene -> mRNA/transcript -> exon/CDS)
  gene_label(g)                     'LOC146345366 (uncharacterized protein ...)' style label
  overlaps(genes, lo, hi)           genes (and the exon/CDS numbers) that overlap lo-hi
"""
import gzip
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
                    description=a.get("description", ""), biotype=a.get("gene_biotype", f[2]),
                    start=s, end=e, strand=f[6], tx=[])
            elif f[2] in TX_TYPES:
                tx[a.get("ID", f"t{len(tx)}")] = dict(id=a.get("ID", ""), parent=a.get("Parent", ""),
                                                     product=a.get("product", ""), start=s, end=e,
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
