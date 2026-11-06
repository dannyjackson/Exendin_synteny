#!/usr/bin/env python3
"""
Rebuild a gene model from the genome when a TOGA2 projection is wrong.

TOGA2 carries the reference transcript's exon structure onto the query, so if the
reference isoform is mis-annotated, the projection is wrong in the query too. In
this project that happened to the HelCha exendin-3 gene: the only isoform TOGA2
had for the reference (LOC_00012325-mRNA-2) uses a 419 bp alternative first exon
whose CDS call is unsupported, so the HelCha model picked up ~110 residues of
non-coding sequence and missed the real first exon 5.9 kb away.

This script takes a trusted query (a published cDNA, or a protein), aligns it to
the genomic region with BLAST, chains the hits into exons, splices them, and
reports the longest ORF. Writes <name>.cds.fa, <name>.aa.fa and <name>.exons.tsv.

Needs: blastn (or tblastn), samtools. Read-only on the genome.

Usage
  # cDNA query (preferred: defines the real start codon)
  python3 2e_rescue_locus.py --genome $HELCHA_FASTA --contig JANEZZ010000386.1 \
      --query AJ580309.fa --name HelCha_A129_EX3 --outdir rescued

  # protein query, when no cDNA exists
  python3 2e_rescue_locus.py --genome $GENOME --contig CTG --query prot.fa --protein \
      --name Some_gene --outdir rescued

Check the printed exon table before trusting the result: the hits should be
collinear, on one strand, and cover the query end to end.
"""
import argparse, os, re, subprocess, sys

B, A = "TCAG", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: A[16 * i + 4 * j + k]
        for i, a in enumerate(B) for j, b in enumerate(B) for k, c in enumerate(B)}


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode:
        sys.exit(f"[error] {cmd[0]} failed:\n{r.stderr.strip()}")
    return r.stdout


def faidx(genome, region, revcomp=False):
    cmd = ["samtools", "faidx"] + (["-i"] if revcomp else []) + [genome, region]
    return "".join(l.strip() for l in run(cmd).splitlines()[1:]).upper()


def translate(nt):
    return "".join(CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def chain(hits, protein):
    """Keep collinear, same-strand hits and trim overlaps in query coordinates."""
    plus = sum(1 for h in hits if h["send"] > h["sstart"])
    strand = "+" if plus >= len(hits) - plus else "-"
    hits = [h for h in hits if (h["send"] > h["sstart"]) == (strand == "+")]
    hits.sort(key=lambda h: (h["qstart"], -h["bits"]))
    step = 3 if protein else 1          # one query residue = 3 genomic bases
    kept, last_q, last_s = [], 0, None
    for h in hits:
        qs, qe, ss, se = h["qstart"], h["qend"], h["sstart"], h["send"]
        if qe <= last_q:
            continue                    # wholly inside a kept hit
        if qs <= last_q:                # trim the 5' overlap at the splice junction
            d = (last_q - qs + 1) * step
            qs = last_q + 1
            ss = ss + d if strand == "+" else ss - d
        if last_s is not None:          # must stay collinear along the genome
            if (strand == "+" and ss <= last_s) or (strand == "-" and ss >= last_s):
                print(f"[warn] skipping non-collinear hit q{qs}-{qe} s{ss}-{se}", file=sys.stderr)
                continue
        kept.append(dict(h, qstart=qs, qend=qe, sstart=ss, send=se))
        last_q, last_s = qe, se
    return strand, kept


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True, help="indexed genome FASTA")
    ap.add_argument("--contig", required=True)
    ap.add_argument("--region", help="start-end within the contig (default: whole contig)")
    ap.add_argument("--query", required=True, help="cDNA FASTA, or protein with --protein")
    ap.add_argument("--protein", action="store_true", help="query is protein; use tblastn")
    ap.add_argument("--name", required=True, help="name for the output sequences")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-identity", type=float, default=80.0)
    ap.add_argument("--min-length", type=int, default=20)
    ap.add_argument("--evalue", default="1e-5")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    region = f"{a.contig}:{a.region}" if a.region else a.contig
    offset = int(a.region.split("-")[0]) - 1 if a.region else 0
    target = os.path.join(a.outdir, f"{a.name}.target.fa")
    with open(target, "w") as fh:
        fh.write(run(["samtools", "faidx", a.genome, region]))

    prog = "tblastn" if a.protein else "blastn"
    cols = "pident length qstart qend sstart send evalue bitscore"
    cmd = [prog, "-query", a.query, "-subject", target, "-outfmt", f"6 {cols}",
           "-evalue", a.evalue]
    # match 2f_scan_interval.py: low-complexity filtering suppresses the short hits
    # that single exons of these small precursors produce
    cmd += ["-seg", "no"] if a.protein else ["-task", "blastn"]
    hits = []
    for line in run(cmd).splitlines():
        f = line.split("\t")
        h = dict(pident=float(f[0]), length=int(f[1]), qstart=int(f[2]), qend=int(f[3]),
                 sstart=int(f[4]), send=int(f[5]), evalue=f[6], bits=float(f[7]))
        if h["pident"] >= a.min_identity and h["length"] >= a.min_length:
            hits.append(h)
    if not hits:
        sys.exit("[error] no BLAST hits passed the filters; widen --region or lower --min-identity")

    strand, exons = chain(hits, a.protein)
    if not exons:
        sys.exit("[error] no collinear chain could be built")

    seq, rows = "", []
    for i, h in enumerate(exons, 1):
        lo, hi = sorted((h["sstart"], h["send"]))
        s = faidx(a.genome, f"{a.contig}:{lo + offset}-{hi + offset}", revcomp=(strand == "-"))
        seq += s
        rows.append([i, h["qstart"], h["qend"], lo + offset, hi + offset, strand,
                     len(s), f"{h['pident']:.1f}"])

    orfs = [m.group(0) for f in range(3)
            for m in re.finditer(r"ATG(?:[ACGT]{3})*?(?:TAA|TAG|TGA)", seq[f:])]
    if not orfs:
        sys.exit(f"[error] spliced {len(seq)} bp but found no complete ORF; check the exon table")
    cds = max(orfs, key=len)
    aa = translate(cds).rstrip("*")

    with open(os.path.join(a.outdir, f"{a.name}.exons.tsv"), "w") as fh:
        fh.write("exon\tq_start\tq_end\tg_start\tg_end\tstrand\tlength\tpident\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(a.outdir, f"{a.name}.cds.fa"), "w") as fh:
        fh.write(f">{a.name}\n{cds[:-3]}\n")          # stop codon removed
    with open(os.path.join(a.outdir, f"{a.name}.aa.fa"), "w") as fh:
        fh.write(f">{a.name}\n{aa}\n")

    print(f"[info] {len(exons)} exons, strand {strand}; spliced {len(seq)} bp; "
          f"CDS {len(cds)} bp; protein {len(aa)} aa", file=sys.stderr)
    print("exon\tq_start\tq_end\tg_start\tg_end\tstrand\tlength\tpident", file=sys.stderr)
    for r in rows:
        print("\t".join(map(str, r)), file=sys.stderr)
    print(f">{a.name}\n{aa}")


if __name__ == "__main__":
    main()
