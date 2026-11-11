#!/usr/bin/env python3
"""
One figure for 8a: each published exendin against every exendin-family gene of HelSus
and HelCha.

  Top    per peptide: where and from which species it was reported, and the best gene of
         each genome with its identity (cell shade = identity; an exact match is marked).
  Below  one block per peptide: the published peptide, then one row per locus (HelSus
         loci, then HelCha loci, each by identity), each cell the residue that locus's
         codon for the peptide residue encodes (tinted = identical, red letter =
         different, red box = stop, dot = not found, | = exon boundary), the genome
         residues before residue 1 and after the last in grey, identity on the right, the
         best locus of each genome starred.

Usage
  python3 8b_plot_published_vs_loci.py --match-dir match --out published_vs_loci
"""
import argparse, csv, importlib.util, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

_sp = importlib.util.spec_from_file_location("plot6b", os.path.join(HERE, "6b_plot_published_map.py"))
P6 = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(P6)
Svg, tint, text_w = P6.Svg, P6.tint, P6.text_w
RED = "#C0392B"
GREY = "#9A9A9A"
SPECIES_NAME = {"HelSus": "H. suspectum", "HelCha": "H. charlesbogerti (H. horridum group)"}


def read_tsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def name_colour(locus):
    """Colour of a step-5 name: EX1-EX4 in the exendin colours, GCGL green, letters grey."""
    import re
    suf = locus.split("_", 1)[1] if "_" in locus else locus
    m = re.match(r"EX(\d)", suf)
    if m:
        return N.FIXED_COLOURS.get(f"exendin-{m.group(1)}", "#555")
    if suf.startswith("GCGL"):
        return N.FIXED_COLOURS.get("GCGL (glucagon-like)", "#009E73")
    return "#777777"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--match-dir", required=True, help="8a output")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default=None)
    ap.add_argument("--min-identity", type=float, default=0.0,
                    help="hide loci below this fraction identity to a peptide (0 = show all) [0]")
    ap.add_argument("--show-no-hit", action="store_true", help="draw a row for loci the peptide does not hit "
                    "(default: listed in one line)")
    ap.add_argument("--width", type=int, default=1320)
    ap.add_argument("--png-width", type=int, default=2640)
    a = ap.parse_args()
    d = a.match_dir
    for f in ("matches.tsv", "residues.tsv", "best.tsv", "peptides.tsv"):
        if not os.path.exists(os.path.join(d, f)):
            sys.exit(f"[error] {d}/{f} missing: run 8a first")
    peps = read_tsv(os.path.join(d, "peptides.tsv"))
    matches = read_tsv(os.path.join(d, "matches.tsv"))
    best = {r["peptide"]: r for r in read_tsv(os.path.join(d, "best.tsv"))}
    res = defaultdict(dict)
    for r in read_tsv(os.path.join(d, "residues.tsv")):
        res[(r["peptide"], r["species"], r["locus"])][int(r["qpos"])] = r
    species = []
    for m in matches:
        if m["species"] not in species:
            species.append(m["species"])

    W, M = a.width, 34
    LX = M + 150 + 30                                # first cell (after locus label and flank)
    nmax = max(len(p["seq"]) for p in peps)
    cw = min(16.0, (W - LX - M - 170 - 4 * 14) / nmax)
    fs = min(12.5, cw * 0.82)
    svg = Svg()
    title = a.title or "Published exendins against every exendin-family gene of the two Heloderma genomes"
    svg.text(M, 32, title, size=18, weight="bold")
    svg.text(M, 52, "Which genome encodes each venom peptide exactly? " +
             " · ".join(f"{sp} = {SPECIES_NAME.get(sp, sp)}" for sp in species), size=11, fill="#555")
    y = 84
    # ---------------- summary ----------------
    cols = [("peptide", 200), ("reported from", 400)] + [(sp, 170) for sp in species] + [("best genome", 120)]
    x = M
    for h, w in cols:
        svg.text(x + 4, y, h, size=11, weight="bold", fill="#333")
        x += w
    svg.line(M, y + 6, M + sum(w for _, w in cols), y + 6, stroke="#999")
    y += 24
    for p in peps:
        b = best.get(p["peptide"], {})
        col = N.FIXED_COLOURS.get(p["group"], "#555")
        x = M
        svg.text(x + 4, y, p["label"], size=11.5, weight="bold", fill=col)
        x += cols[0][1]
        import re
        src = p["source"]
        who = re.match(r"([A-Z][a-z]+)( et al\.)? (\d{4})", src)
        frm = re.search(r"from (H\. [a-z]+) venom", src)
        txt = (f"{who.group(1)} {who.group(3)}: " if who else "") + (f"{frm.group(1)} venom" if frm else src[:40])
        if "typical of H. horridum" in src:
            txt += "; typical of H. horridum (Vandermeers 1987)"
        svg.text(x + 4, y, txt, size=10, fill="#444")
        x += cols[1][1]
        L = int(p["seq"] and len(p["seq"]))
        for sp in species:
            k = int(b.get(f"{sp}_identical", 0) or 0)
            f = k / L if L else 0
            svg.rect(x + 2, y - 13, cols[2][1] - 8, 18, tint(col, 0.12 + 0.75 * f ** 3), rx=3)
            lab = f"{k}/{L} ({100 * f:.0f}%)  {b.get(f'{sp}_locus', '-')}"
            svg.text(x + 8, y, ("✓ " if k == L else "") + lab, size=10, weight="bold" if k == L else "normal",
                     fill="#1E1A1D")
            x += 170
        bg = b.get("best_genome", "")
        mg = b.get("margin", "")
        svg.text(x + 4, y, f"{bg} (+{mg})" if bg and mg not in ("", None) else bg or "-", size=11, weight="bold")
        y += 24
    y += 14
    svg.text(M, y, "Shade = identity of the best gene in each genome (exact match ✓). +n = how many more residues "
                   "the best genome's gene matches than the other genome's best.", size=9.5, fill="#666", italic=True)
    y += 30

    # ---------------- per-peptide blocks ----------------
    for p in peps:
        pid, seq = p["peptide"], p["seq"]
        col = N.FIXED_COLOURS.get(p["group"], "#555")
        Lp = len(seq)
        svg.rect(M, y - 16, W - 2 * M, 44, tint(col, 0.15), rx=4)
        svg.text(M + 8, y, p["label"], size=13, weight="bold", fill=col)
        y += 20
        svg.text(LX - 6, y, "published", size=10, anchor="end", fill="#444")
        for q in range(1, Lp + 1):
            svg.text(LX + (q - 1) * cw + cw / 2, y, seq[q - 1], size=fs, anchor="middle", mono=True, weight="bold")
        y += 22
        for sp in species:
            ms = [m for m in matches if m["peptide"] == pid and m["species"] == sp
                  and int(m["identical"]) / Lp >= a.min_identity and (a.show_no_hit or m["note"] == "")]
            nohit = [m["locus"] for m in matches if m["peptide"] == pid and m["species"] == sp and m["note"] != ""]
            ms.sort(key=lambda m: (-int(m["identical"]), m["locus"]))
            svg.text(M, y + 2, SPECIES_NAME.get(sp, sp), size=9.5, fill="#666", italic=True)
            y += 15
            for m in ms:
                rr = res.get((pid, sp, m["locus"]), {})
                star = m.get("best_in_species") == "yes"
                lc = name_colour(m["locus"])
                svg.text(M + 10, y, ("★ " if star else "  ") + m["locus"], size=11,
                         weight="bold" if star else "normal", fill=lc)
                bf = "" if m["before"] in ("-", "") else m["before"][-2:]
                for k, ch in enumerate(bf):
                    svg.text(LX - (len(bf) - k) * cw + cw / 2 - 2, y, ch, size=fs * 0.9, anchor="middle",
                             mono=True, fill=GREY)
                prev_ex = None
                last_x = LX
                for q in range(1, Lp + 1):
                    x = LX + (q - 1) * cw
                    r = rr.get(q)
                    if r and r["exon"] and prev_ex and r["exon"] != prev_ex:
                        svg.line(x, y - cw + 1, x, y + 4, stroke="#333", sw=1.4)
                    if r and r["exon"]:
                        prev_ex = r["exon"]
                    if not r or r["status"] == "unmapped" or r["codon_start"] == "":
                        svg.text(x + cw / 2, y, "·", size=fs, anchor="middle", mono=True, fill="#BBB")
                    elif r["status"] == "identical":
                        svg.rect(x + 0.5, y - cw + 3, cw - 1, cw, tint(col, 0.3))
                        svg.text(x + cw / 2, y, r["gres"], size=fs, anchor="middle", mono=True)
                    elif r["status"] == "stop":
                        svg.rect(x + 0.5, y - cw + 3, cw - 1, cw, RED)
                        svg.text(x + cw / 2, y, "*", size=fs, anchor="middle", mono=True, fill="white")
                    elif r["status"] == "gap":
                        svg.text(x + cw / 2, y, "–", size=fs, anchor="middle", mono=True, fill=RED)
                    else:
                        svg.text(x + cw / 2, y, r["gres"], size=fs, anchor="middle", mono=True, fill=RED,
                                 weight="bold")
                    last_x = x + cw
                af = "" if m["after"] in ("-", "") else m["after"]
                for k, ch in enumerate(af):
                    svg.text(last_x + 2 + k * cw + cw / 2, y, ch, size=fs * 0.9, anchor="middle", mono=True, fill=GREY)
                k_, f_ = int(m["identical"]), int(m["identical"]) / Lp
                svg.text(W - M, y, f"{k_}/{Lp}  {100 * f_:.0f}%" if m["note"] == "" else "no hit", size=11,
                         anchor="end", weight="bold" if star else "normal", fill="#1E1A1D" if f_ >= 0.8 else "#777")
                y += cw + 6
            if nohit and not a.show_no_hit:
                svg.text(M + 10, y - 2, "no hit: " + ", ".join(nohit), size=9.5, fill="#999", italic=True)
                y += 14
            y += 4
        y += 16
    # legend
    lx = M
    for fill, lab in [(tint("#888888", 0.3), "identical"), ("#FFFFFF", "different (red letter)"), (RED, "stop codon"),
                      ("#EEEEEE", "· not found"), ("#333333", "| exon boundary"),
                      ("#FFFFFF", "grey letters: genome before residue 1 / after the last"),
                      ("#FFFFFF", "★ best gene of that genome")]:
        svg.rect(lx, y - 9, 10, 10, fill, stroke="#999", sw=0.6)
        svg.text(lx + 14, y, lab, size=9.5, fill="#444")
        lx += 30 + text_w(lab, 9.5)
        if lx > W - 250:
            lx = M; y += 16
    y += 16
    svg.text(M, y, "Locus names: step 5 (EX1–EX4 = published-peptide match, GCGL, letters = position, coloured "
                   "accordingly). Mapping: 8a (tblastn per locus, splice-aware chaining as 7a).", size=9.5, fill="#777")
    H = int(y + 20)
    out = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">\n'
           f'<rect width="{W}" height="{H}" fill="white"/>\n' + "\n".join(svg.el) + "\n</svg>\n")
    open(a.out + ".svg", "w").write(out)
    N.load_engine().write_png(out, a.out + ".png", a.png_width)
    N.info(f"wrote {a.out}.svg")


if __name__ == "__main__":
    main()
