#!/usr/bin/env python3
"""
CDS of outgroup genes (GCG, GIP, GCGL ...) for the exendin tree, from two kinds of source,
written as one FASTA to add with 10b_apply_curation.py --add-cds (10_toga2_exendin_tree.md).

  --gene CODE:SYMBOL[,SYMBOL]   a gene from that species' NCBI annotation in ${GENOMEDIR}
                                (species.tsv + 1b layout: genomic.gff + *_genomic.fna).
                                The gene is found by symbol (gene=/Name=), or else by
                                --describe (description regex). Its transcript tagged
                                'MANE Select' is used if there is one, else the transcript
                                with the longest CDS. The CDS is cut from the genome
                                (samtools faidx) and reverse-complemented on the - strand.
  --cds-fasta CODE=FILE[:NAME]  CDS records already fetched, e.g. by
                                efetch -db nuccore -id NM_002054 -format fasta_cds_na.
                                The tip is CODE_<gene> from the record's [gene=...] tag, or
                                CODE_NAME if NAME is given (for records without a gene tag).

Every CDS is checked: length a multiple of 3, starts ATG, internal stops, final stop
(removed). Problems are reported, not fixed.

Tip names: <CODE>_<SYMBOL>, e.g. GalGal_GCG, HomSap_GIP, GalGal_GCGL.

Writes --out (CDS FASTA) and --out.tsv (source, transcript, length, checks, protein).

Usage
  python3 10d_outgroup_cds.py --species-table ../species.tsv --genomedir ${GENOMEDIR} \
      --gene GalGal:GCG,GIP --gene HelSus:GCG,GIP \
      --cds-fasta HomSap=outgroups/human_cds.fa --cds-fasta GalGal=outgroups/EU718628_cds.fa:GCGL \
      --out outgroups/outgroups_cds.fa
"""
import argparse, gzip, os, re, subprocess, sys
from collections import defaultdict
from urllib.parse import unquote

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

_B, _AA = "TCAG", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: _AA[16 * i + 4 * j + k] for i, a in enumerate(_B) for j, b in enumerate(_B) for k, c in enumerate(_B)}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
DESCRIBE = {"GCG": r"^glucagon$|proglucagon", "GIP": r"gastric inhibitory polypeptide|glucose-dependent insulinotropic polypeptide$",
            "GCGL": r"glucagon[- ]like(?! peptide [12] receptor)"}


def translate(nt):
    return "".join(CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def gff_models(gff):
    genes, tx, cds = {}, {}, defaultdict(list)
    with opener(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 9:
                continue
            at = {k: unquote(v) for k, v in (kv.split("=", 1) for kv in t[8].strip().strip(";").split(";") if "=" in kv)}
            if t[2] in ("gene", "pseudogene"):
                genes[at.get("ID", "")] = dict(symbol=at.get("gene") or at.get("Name", ""),
                                               desc=at.get("description", ""), biotype=at.get("gene_biotype", ""))
            elif t[2] in ("mRNA", "transcript"):
                tx[at.get("ID", "")] = dict(parent=at.get("Parent", "").split(",")[0], name=at.get("transcript_id") or at.get("Name") or at.get("ID"),
                                            tag=at.get("tag", ""), product=at.get("product", ""))
            elif t[2] == "CDS":
                cds[at.get("Parent", "").split(",")[0]].append((t[0], int(t[3]), int(t[4]), t[6], int(t[7]) if t[7].isdigit() else 0))
    return genes, tx, cds


def extract(fasta, parts):
    parts = sorted(parts, key=lambda p: p[1])
    strand = parts[0][3]
    seq = ""
    for c, s, e, *_ in parts:
        out = subprocess.run(["samtools", "faidx", fasta, f"{c}:{s}-{e}"], capture_output=True, text=True)
        if out.returncode:
            sys.exit(f"[error] samtools faidx {c}:{s}-{e}: {out.stderr.strip()}")
        seq += "".join(out.stdout.splitlines()[1:]).upper()
    if strand == "-":
        seq = seq.translate(COMP)[::-1]
    return seq


def read_cds_fasta(path):
    recs, cur = [], None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            cur = [line[1:], []]; recs.append(cur)
        elif cur:
            cur[1].append(line.strip().upper())
    return [(h, "".join(s)) for h, s in recs]


def check(nt):
    notes = []
    if len(nt) % 3:
        notes.append("length not a multiple of 3")
    if not nt.startswith("ATG"):
        notes.append("no ATG start (partial CDS?)")
    aa = translate(nt)
    if aa.endswith("*"):
        nt, aa = nt[:-3], aa[:-1]
    else:
        notes.append("no final stop (partial CDS?)")
    if "*" in aa:
        notes.append(f"{aa.count('*')} internal stop(s)")
    return nt, aa, notes


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--species-table", help="species.tsv (code -> accession) for --gene")
    ap.add_argument("--genomedir", help="${GENOMEDIR}")
    ap.add_argument("--gene", action="append", default=[], help="CODE:SYMBOL[,SYMBOL] from that species' NCBI annotation")
    ap.add_argument("--describe", action="append", default=[], help="SYMBOL=REGEX: description to match if the symbol is absent")
    ap.add_argument("--cds-fasta", action="append", default=[], help="CODE=FILE[:NAME] of fetched CDS records")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    describe = dict(DESCRIBE)
    for d in a.describe:
        k, _, v = d.partition("=")
        describe[k] = v
    rows = []

    if a.gene:
        if not (a.species_table and a.genomedir):
            sys.exit("[error] --gene needs --species-table and --genomedir")
        acc = {r["code"]: r["accession"] for r in N.read_species_table(a.species_table)}
        for spec in a.gene:
            code, _, syms = spec.partition(":")
            if code not in acc or acc[code] in ("", "NA"):
                sys.exit(f"[error] {code}: no accession in {a.species_table}")
            f = N.find_genome_files(a.genomedir, code, acc[code])
            if not (f["gff"] and f["fasta"]):
                sys.exit(f"[error] {code}: GFF or genome FASTA missing under {a.genomedir} (run 1b)")
            if not os.path.exists(f["fasta"] + ".fai"):
                subprocess.run(["samtools", "faidx", f["fasta"]], check=True)
            genes, tx, cds = gff_models(f["gff"])
            for sym in syms.split(","):
                hits = [gid for gid, g in genes.items() if g["symbol"].upper() == sym.upper()]
                how = "symbol"
                if not hits and sym.upper() in describe:
                    rx = re.compile(describe[sym.upper()], re.I)
                    hits = [gid for gid, g in genes.items() if rx.search(g["desc"]) and g["biotype"] in ("protein_coding", "")]
                    how = f"description /{describe[sym.upper()]}/"
                if not hits:
                    N.warn(f"{code}: no gene {sym} (by symbol or description); skipped")
                    continue
                if len(hits) > 1:
                    N.warn(f"{code} {sym}: {len(hits)} genes match by {how}: "
                           + ", ".join(f"{genes[h]['symbol']} ({genes[h]['desc'][:40]})" for h in hits[:5]) + "; using the first")
                gid = hits[0]
                txs = [(tid, t) for tid, t in tx.items() if t["parent"] == gid and cds.get(tid)]
                if not txs:
                    N.warn(f"{code} {sym}: gene {gid} has no coding transcript; skipped")
                    continue
                mane = [x for x in txs if "MANE" in x[1]["tag"]]
                tid, t = (mane or sorted(txs, key=lambda x: -sum(e - s + 1 for _, s, e, *_ in cds[x[0]])))[0]
                nt, aa, notes = check(extract(f["fasta"], cds[tid]))
                rows.append([f"{code}_{sym.upper()}", f"{code} {acc[code]} gene {genes[gid]['symbol']} ({how})",
                             t["name"], nt, aa, notes, len(txs)])

    for spec in a.cds_fasta:
        code, _, rest = spec.partition("=")
        path, _, name = rest.partition(":")
        for head, nt in read_cds_fasta(path):
            g = re.search(r"\[gene=([^\]]+)\]", head)
            sym = name or (g.group(1) if g else None)
            if not sym:
                N.warn(f"{path}: record '{head[:60]}' has no [gene=] tag; give a name as {code}={path}:NAME")
                continue
            tid = re.search(r"\[protein_id=([^\]]+)\]", head)
            nt, aa, notes = check(nt)
            rows.append([f"{code}_{sym.upper()}", f"{os.path.basename(path)}", head.split()[0] + (f" ({tid.group(1)})" if tid else ""),
                         nt, aa, notes, 1])

    seen = defaultdict(int)
    with open(a.out, "w") as fo, open(a.out + ".tsv", "w") as ft:
        ft.write("tip\tsource\ttranscript\tcds_bp\taa_len\tchecks\tn_coding_transcripts\tprotein\n")
        for tip, src, tr, nt, aa, notes, ntx in rows:
            seen[tip] += 1
            if seen[tip] > 1:
                N.warn(f"{tip}: more than one record; the extra one is named {tip}_{seen[tip]}")
                tip = f"{tip}_{seen[tip]}"
            fo.write(f">{tip}\n{nt}\n")
            ft.write(f"{tip}\t{src}\t{tr}\t{len(nt)}\t{len(aa)}\t{'; '.join(notes) or 'ok'}\t{ntx}\t{aa}\n")
            N.info(f"{tip}: {tr}, {len(aa)} aa{'; ' + '; '.join(notes) if notes else ''}  {aa[:40]}...")
    N.info(f"wrote {len(rows)} CDS to {a.out}")


if __name__ == "__main__":
    main()
