#!/usr/bin/env python3
"""
Apply the curation table to the collected exendin CDS, so every curation decision is a line
in a file rather than a hand edit to a FASTA (10_toga2_exendin_tree.md).

curation.tsv (written by 10a with every tip on 'keep'; '#' lines are comments):
  tip  action  value  note
  keep     use as is
  drop     leave out
  rename   value = new tip name
  trim     value = first-last amino acid to KEEP, 1-based in the unaligned protein (the
           CDS is cut to match), e.g. 43-129
  Several rows for one tip are applied in order (trim, then rename).
A tip with no row is kept and reported. --add-cds adds sequences from outside TOGA2
(e.g. a rebuilt model); they are logged as 'added'.

Writes --out-prefix _cds.fa, _aa.fa (internal stops as X) and _log.tsv (what was done to
each tip, and the note).

Usage
  python3 10b_apply_curation.py --cds collect/exendin_cds.fa --curation curation.tsv \
      --out-prefix curated/exendin 2> curated.log
"""
import argparse, importlib.util, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("toga2_collect", os.path.join(HERE, "3b_collect_toga2_exendins.py"))
T3B = importlib.util.module_from_spec(spec); spec.loader.exec_module(T3B)


def read_fa(path):
    out, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; out[name] = []
        elif name:
            out[name].append(line.upper())
    return {k: "".join(v) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cds", required=True, help="exendin_cds.fa from 10a")
    ap.add_argument("--curation", required=True)
    ap.add_argument("--add-cds", action="append", default=[], help="FASTA of extra CDS to add")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    seqs = read_fa(a.cds)
    origin = {k: "TOGA2" for k in seqs}
    for f in a.add_cds:
        for k, v in read_fa(f).items():
            if k in seqs:
                sys.exit(f"[error] {k} from {f} is already a tip")
            seqs[k], origin[k] = v, f"added:{os.path.basename(f)}"
    rules = defaultdict(list)
    for line in open(a.curation):
        if line.startswith("#") or not line.strip():
            continue
        t = (line.rstrip("\n").split("\t") + ["", "", ""])[:4]
        if t[0] == "tip":
            continue
        if t[0] not in seqs:
            print(f"[warn] curation row for {t[0]}: not a tip; ignored", file=sys.stderr)
            continue
        rules[t[0]].append((t[1].strip().lower(), t[2].strip(), t[3].strip()))

    out, log, names = [], [], set()
    for tip, nt in seqs.items():
        name, done, notes, drop = tip, [], [], False
        if tip not in rules and origin[tip] == "TOGA2":
            print(f"[warn] {tip}: no curation row; kept", file=sys.stderr)
        for act, val, note in rules.get(tip, [("keep", "", "")]):
            if note:
                notes.append(note)
            if act == "keep":
                done.append("keep")
            elif act == "drop":
                drop = True; done.append("drop")
            elif act == "rename":
                name = val; done.append(f"rename->{val}")
            elif act == "trim":
                s, _, e = val.partition("-")
                s, e = int(s), int(e)
                if not 1 <= s <= e <= len(nt) // 3:
                    sys.exit(f"[error] {tip}: trim {val} outside 1-{len(nt) // 3}")
                nt = nt[(s - 1) * 3:e * 3]; done.append(f"trim {val}")
            else:
                sys.exit(f"[error] {tip}: unknown action '{act}'")
        if not drop:
            if name in names:
                sys.exit(f"[error] two tips would be called {name}")
            names.add(name)
            out.append((name, nt))
        aa = T3B.translate(nt)
        log.append([tip, name if not drop else "-", origin[tip], ";".join(done), len(nt), len(aa),
                    aa.count("*"), " | ".join(notes)])

    with open(a.out_prefix + "_cds.fa", "w") as fn, open(a.out_prefix + "_aa.fa", "w") as fa:
        for n, nt in out:
            fn.write(f">{n}\n{nt}\n")
            fa.write(f">{n}\n{T3B.translate(nt).replace('*', 'X')}\n")
    with open(a.out_prefix + "_log.tsv", "w") as fh:
        fh.write("tip\tnew_name\torigin\tactions\tcds_bp\taa_len\tinternal_stops\tnote\n")
        for r in log:
            fh.write("\t".join(map(str, r)) + "\n")
    nd = sum(";drop" in ";" + r[3] for r in log)
    print(f"[info] {len(seqs)} tips in, {nd} dropped, {len(out)} written to {a.out_prefix}_cds.fa / _aa.fa", file=sys.stderr)
    for r in log:
        if r[3] not in ("keep", ""):
            print(f"[info]   {r[0]}: {r[3]}" + (f"  ({r[7]})" if r[7] else ""), file=sys.stderr)


if __name__ == "__main__":
    main()
