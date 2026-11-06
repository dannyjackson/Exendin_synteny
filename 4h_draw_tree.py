#!/usr/bin/env python3
"""
Draw the step-4 gene tree with the same colour per name group as the synteny plots.

This runs lib_draw_tree.py (the old repo's 2c_draw_tree.py, copied unchanged), after
replacing its fixed colour table with lib_exendin_names.FIXED_COLOURS, so exendin-2 is
the same colour in the tree and in all three synteny figures. All options are 2c's:
  --tree --groups --outgroup --title --min-support --width --row --png --png-scale
Node labels are SH-aLRT/UFBoot; nodes with UFBoot >= --min-support are dark. Bold =
Heloderma.

Usage
  python3 4h_draw_tree.py --tree tree/exendin_aa.rooted.nwk --groups tree/groups.tsv \
      --outgroup $(cat tree/outgroup.txt) --title "..." --png tree/exendin_aa.png > tree/exendin_aa.svg
"""
import importlib.util, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

spec = importlib.util.spec_from_file_location("draw_tree", os.path.join(HERE, "lib_draw_tree.py"))
T = importlib.util.module_from_spec(spec)
spec.loader.exec_module(T)
T.FIXED = dict(N.FIXED_COLOURS, outgroup="#333333")
T.PALETTE = N.EXTRA_COLOURS

if __name__ == "__main__":
    T.main()
