#!/usr/bin/env python3
"""
Figure for 7a: is each locus an intact copy of the helodermin (exendin-2) precursor?

  A  One row per locus, in precursor coordinates (residue 1 to the end of P04204), under
     the precursor itself (signal peptide, propeptide, mature helodermin marked). Each
     cell is the residue the locus's codon for that precursor residue encodes:
     tinted = identical, red letter = different, red box = stop codon, dot = not found.
     Between cells: exon boundaries (intron length, donor..acceptor; red if no splice
     site keeps the frame), frameshifts (purple, with the indel), and the end of the
     assembled ORF (black bar). The verdict of 7a is on the left.
  B  The 7a summary: verdict and reasons, exons, splice sites, start codon, identity to
     the precursor overall and per region, the KR site before the mature peptide, TOGA2.
  C  dN/dS of each --pair (Nei-Gojobori over codons placed on the same precursor residue).

Usage
  python3 7b_plot_pseudogene_check.py --check-dir check --out pseudogene_check
Writes <out>.svg and <out>.png.
"""
import argparse, csv, importlib.util, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

_spec = importlib.util.spec_from_file_location("plot6b", os.path.join(HERE, "6b_plot_published_map.py"))
P6 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(P6)
Svg, tint, text_w = P6.Svg, P6.tint, P6.text_w

VERDICT = {"intact": "#2E8B57", "disrupted": "#C0392B", "partial": "#D4A017", "not found": "#888888"}
FS = "#7B3FA0"
RED = "#C0392B"


def read_tsv(p):
    if not os.path.exists(p):
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check-dir", required=True, help="7a output")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default=None)
    ap.add_argument("--width", type=int, default=1700)
    ap.add_argument("--png-width", type=int, default=3400)
    a = ap.parse_args()
    d = a.check_dir
    summ = read_tsv(os.path.join(d, "summary.tsv"))
    if not summ:
        sys.exit(f"[error] no {d}/summary.tsv: run 7a first")
    exons = defaultdict(list)
    for r in read_tsv(os.path.join(d, "exons.tsv")):
        exons[r["locus"]].append(r)
    events = defaultdict(list)
    for r in read_tsv(os.path.join(d, "events.tsv")):
        events[r["locus"]].append(r)
    res = defaultdict(dict)
    for r in read_tsv(os.path.join(d, "residues.tsv")):
        res[r["locus"]][int(r["qpos"])] = r
    regions = {r["region"]: (int(r["start"]), int(r["end"]), r["sequence"]) for r in read_tsv(os.path.join(d, "reference.tsv"))}
    pairs = read_tsv(os.path.join(d, "pairs.tsv"))
    prec = "".join(regions[k][2] for k in ("signal", "propeptide", "mature") if k in regions)
    Lq = len(prec)
    # rows: by species, control before test
    order = sorted(summ, key=lambda r: (r.get("species", ""), r.get("role") != "control", r["name"]))

    W, M = a.width, 36
    LX = M + 250                                   # first residue cell
    cw = min(14.0, (W - LX - M - 10) / Lq)
    fs = min(12.0, cw * 0.88)

    def X(q):                                      # left edge of residue q's cell
        return LX + (q - 1) * cw
    svg = Svg()
    title = a.title or "Exendin-2-like loci against the helodermin precursor (P04204): intact gene or pseudogene?"
    svg.text(M, 32, title, size=19, weight="bold")
    y = 70
    svg.text(M, y, "A", size=22, weight="bold")
    svg.text(M + 26, y - 1, "Each locus, codon by codon, in precursor coordinates", size=14, weight="bold")
    y += 28
    # precursor regions
    rc = {"signal": "#B8B0B6", "propeptide": "#E3D9C6", "mature": N.FIXED_COLOURS.get("exendin-2", "#D55E00")}
    for k, (r0, r1, _) in regions.items():
        svg.rect(X(r0), y, (r1 - r0 + 1) * cw, 12, rc.get(k, "#DDD"), rx=2)
        lab = {"signal": "signal peptide", "propeptide": "propeptide", "mature": "mature helodermin (exendin-2)"}.get(k, k)
        svg.text(X(r0) + (r1 - r0 + 1) * cw / 2, y - 4, f"{lab} {r0}–{r1}", size=10, anchor="middle", fill="#444")
    y += 30
    svg.text(LX - 10, y, "precursor", size=11, anchor="end", weight="bold")
    for q in range(1, Lq + 1):
        svg.text(X(q) + cw / 2, y, prec[q - 1], size=fs, anchor="middle", mono=True, weight="bold")
        if q % 10 == 0:
            svg.text(X(q) + cw / 2, y - 14, str(q), size=8, anchor="middle", fill="#999")
    y += 16
    ex2c = N.FIXED_COLOURS.get("exendin-2", "#D55E00")
    for s in order:
        nm = s["name"]
        y += 46
        v = s.get("verdict", "?")
        vc = VERDICT.get(v, "#888")
        svg.text(M, y, nm, size=13, weight="bold")
        svg.text(M, y + 15, f"{s.get('role', '')} · {s.get('contig', '')}", size=9.5, fill="#666")
        svg.rect(M + 150, y - 12, 82, 18, vc, rx=9)
        svg.text(M + 191, y + 1, v, size=11, anchor="middle", weight="bold", fill="white")
        # residue cells
        rr = res.get(nm, {})
        for q in range(1, Lq + 1):
            r = rr.get(q)
            x = X(q)
            if not r or r["status"] == "unmapped" or r["codon_start"] == "":
                svg.text(x + cw / 2, y, "·", size=fs, anchor="middle", mono=True, fill="#BBB")
                continue
            st, g = r["status"], r["gres"]
            if st == "identical":
                svg.rect(x + 0.5, y - cw + 2, cw - 1, cw + 1, tint(ex2c, 0.28))
                svg.text(x + cw / 2, y, g, size=fs, anchor="middle", mono=True)
            elif st == "stop":
                svg.rect(x + 0.5, y - cw + 2, cw - 1, cw + 1, RED)
                svg.text(x + cw / 2, y, "*", size=fs, anchor="middle", mono=True, fill="white", weight="bold")
            elif st == "gap":
                svg.text(x + cw / 2, y, "–", size=fs, anchor="middle", mono=True, fill=RED)
            else:
                svg.text(x + cw / 2, y, g, size=fs, anchor="middle", mono=True, fill=RED, weight="bold")
            if r.get("ext") == "yes":
                svg.line(x + 1, y + 3, x + cw - 1, y + 3, stroke="#999", sw=0.8)
        # exon boundaries
        ex = sorted(exons.get(nm, []), key=lambda e: int(e["exon"]))
        for e in ex[:-1]:
            qb = int(e["q_end"])
            xb = X(qb + 1)
            ok = e.get("splice_ok") == "yes"
            col = "#555" if ok else RED
            svg.line(xb, y - cw - 4, xb, y + 6, stroke=col, sw=1.4, dash=None if ok else "3,2")
            il = e.get("intron_after", "")
            lab = (f"{int(il):,} bp " if il.isdigit() else "") + f"{e.get('donor_after', '')}..{e.get('acceptor_after', '')}"
            svg.text(xb, y - cw - 7, lab + ("" if ok else " ✗"), size=8.5, anchor="middle", fill=col,
                     weight="normal" if ok else "bold")
        if not ex:
            svg.text(LX, y, "precursor not found at this locus", size=11, fill="#999", italic=True)
        # events
        for ev in events.get(nm, []):
            q = int(ev["qpos"] or 0)
            if ev["kind"] == "frameshift":
                xb = X(q + 1)
                svg.poly([(xb, y + 6), (xb - 5, y + 15), (xb + 5, y + 15)], fill=FS)
                det = ev["detail"].split(" bp")[0] + " bp" if " bp" in ev["detail"] else "frameshift"
                svg.text(xb, y + 26, det.replace("CDS length", "len"), size=8.5, anchor="middle", fill=FS, weight="bold")
            elif ev["kind"] == "orf_end" and q:
                xb = X(q + 1)
                svg.rect(xb - 1.5, y - cw - 2, 3, cw + 10, "#1E1A1D")
                svg.text(xb + 4, y + 16, "ORF ends", size=8.5, fill="#1E1A1D")
            elif ev["kind"] == "start":
                svg.text(X(1), y + 16, ev["detail"].split("(")[-1].rstrip(")"), size=8.5, fill=RED)
        y += 8
    y += 34
    # legend
    lx = LX
    for fill, lab in [(tint(ex2c, 0.28), "identical to the precursor"), ("#FFFFFF", "different (red letter)"),
                      (RED, "stop codon"), ("#EEEEEE", "· not found"), ("#555", "exon boundary: intron, donor..acceptor"),
                      (RED, "✗ no splice site keeps the frame"), (FS, "▲ frameshift (indel, bp)"),
                      ("#1E1A1D", "end of the assembled ORF")]:
        svg.rect(lx, y - 9, 10, 10, fill, stroke="#999", sw=0.6)
        svg.text(lx + 14, y, lab, size=9.5, fill="#444")
        lx += 30 + text_w(lab, 9.5)
        if lx > W - 260:
            lx = LX; y += 16
    y += 36

    # ================= B: summary table =================
    svg.text(M, y, "B", size=22, weight="bold")
    svg.text(M + 26, y - 1, "Verdict and evidence (7a summary.tsv)", size=14, weight="bold")
    y += 24
    cols = [("locus", "name", 120), ("verdict", "verdict", 80), ("exons", "exons", 48), ("introns (bp)", "introns", 100),
            ("splice", "splice", 110), ("start", "start_codon", 46), ("identity", "identity", 64),
            ("signal", "signal_id", 52), ("propep.", "propeptide_id", 56), ("mature", "mature_id", 56),
            ("KR site", "kr_site", 54), ("TOGA2", "toga2", 150), ("reasons", "reasons", 0)]
    used = sum(c[2] for c in cols)
    cols[-1] = ("reasons", "reasons", W - 2 * M - used)
    x = M
    for h, _, w in cols:
        svg.text(x + 4, y, h, size=10.5, weight="bold", fill="#333")
        x += w
    svg.line(M, y + 5, W - M, y + 5, stroke="#999")
    y += 20
    for s in order:
        x = M
        lines_used = 1
        for h, k, w in cols:
            val = s.get(k, "") or "-"
            col = VERDICT.get(val, "#222") if k == "verdict" else "#222"
            if k == "reasons":
                lim = max(10, int(w / (0.55 * 10)))
                parts, cur = [], ""
                for word in val.split(" "):
                    if len(cur) + len(word) + 1 > lim:
                        parts.append(cur); cur = word
                    else:
                        cur = (cur + " " + word).strip()
                parts.append(cur)
                for i, p in enumerate(parts):
                    svg.text(x + 4, y + 13 * i, p, size=10, fill="#333")
                lines_used = max(lines_used, len(parts))
            else:
                lim = max(4, int(w / (0.55 * 10)))
                svg.text(x + 4, y, val if len(val) <= lim else val[:lim - 1] + "…", size=10, fill=col,
                         weight="bold" if k in ("name", "verdict") else "normal")
            x += w
        y += 13 * lines_used + 7
        svg.line(M, y - 9, W - M, y - 9, stroke="#EEE")
    y += 20

    # ================= C: dN/dS =================
    if pairs:
        svg.text(M, y, "C", size=22, weight="bold")
        svg.text(M + 26, y - 1, "Selection: dN/dS between loci (Nei-Gojobori, codons placed on the same precursor residue)",
                 size=14, weight="bold")
        y += 24
        pc = [("pair", 260), ("codons", 70), ("syn. diffs", 80), ("nonsyn. diffs", 95), ("dS", 70), ("dN", 70),
              ("dN/dS", 70)]
        x = M
        for h, w in pc:
            svg.text(x + 4, y, h, size=10.5, weight="bold", fill="#333")
            x += w
        svg.line(M, y + 5, M + sum(w for _, w in pc), y + 5, stroke="#999")
        y += 20
        for p in pairs:
            def f(k, nd=3):
                try:
                    return f"{float(p[k]):.{nd}f}"
                except (KeyError, ValueError):
                    return "-"
            vals = [f"{p['a']}  vs  {p['b']}", p.get("codons", ""), f("Sd", 1), f("Nd", 1), f("dS"), f("dN"), f("w", 2)]
            x = M
            for (h, w), v in zip(pc, vals):
                bold = h == "dN/dS"
                svg.text(x + 4, y, v, size=10.5, fill="#222", weight="bold" if bold else "normal")
                x += w
            y += 17
        y += 6
        svg.text(M, y, "dN/dS near 1: no constraint on the protein (a pseudogene evolving neutrally). Well below 1: "
                       "purifying selection (a functional gene). With few differences (see syn./nonsyn. diffs) the "
                       "estimate is rough.", size=10, fill="#666", italic=True)
        y += 20
    H = int(y + 20)
    out_svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">\n'
               f'<rect width="{W}" height="{H}" fill="white"/>\n' + "\n".join(svg.el) + "\n</svg>\n")
    open(a.out + ".svg", "w").write(out_svg)
    N.load_engine().write_png(out_svg, a.out + ".png", a.png_width)
    N.info(f"wrote {a.out}.svg")


if __name__ == "__main__":
    main()
