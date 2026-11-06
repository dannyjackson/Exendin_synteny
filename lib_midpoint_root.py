#!/usr/bin/env python3
"""
Midpoint-root an IQ-TREE .treefile, keeping each support value on the branch it
belongs to. No dependencies beyond Python 3.

IQ-TREE writes support (e.g. "85.3/97") as an internal-node label, but it describes
the branch above that node. Naive re-rooting leaves labels on nodes and so moves
them to the wrong branches; here labels are stored on edges and re-attached to the
child side of each edge in the new orientation.

Used by 4_toga2_exendin_tree.md when no outgroup is available from the TOGA2 gene set.

Usage
  python3 4c_midpoint_root.py exendin.codon.treefile > exendin.codon.midpoint.nwk
"""
import re, sys
from collections import defaultdict


def parse(s):
    """Newick -> (adjacency {node: {nbr: [length, label]}}, tip names {node: name})."""
    s = re.sub(r"\[[^\]]*\]", "", s.strip()).rstrip(";")
    adj, names = defaultdict(dict), {}
    nid, stack, cur, i = [0], [], 0, 0
    parent = {0: None}
    def new(par):
        nid[0] += 1
        parent[nid[0]] = par
        return nid[0]
    lab, ln = {}, {}
    while i < len(s):
        c = s[i]
        if c == "(":
            stack.append(cur); cur = new(cur); i += 1
        elif c == ",":
            cur = new(stack[-1]); i += 1
        elif c == ")":
            cur = stack.pop(); i += 1
        elif c == ":":
            m = re.match(r":([^,():;]+)", s[i:]); ln[cur] = float(m.group(1)); i += len(m.group(0))
        else:
            m = re.match(r"('[^']*'|[^,():;]+)", s[i:]); lab[cur] = m.group(1).strip("'"); i += len(m.group(0))
    # node 0 is the root clade; the first "(" opened its first child
    kids = defaultdict(list)
    for n, p in parent.items():
        if p is not None:
            kids[p].append(n)
    root = 0
    for n, p in parent.items():
        if p is None:
            continue
        is_tip = not kids[n]
        edge_lab = "" if is_tip else lab.get(n, "")
        adj[p][n] = adj[n][p] = [ln.get(n, 0.0), edge_lab]
        if is_tip:
            names[n] = lab.get(n, f"tip{n}")
    return adj, names, root


def far(adj, start):
    dist, prev, stack = {start: 0.0}, {start: None}, [start]
    while stack:
        n = stack.pop()
        for m, (l, _) in adj[n].items():
            if m not in dist:
                dist[m], prev[m] = dist[n] + l, n
                stack.append(m)
    end = max(dist, key=dist.get)
    return end, dist, prev


def main():
    adj, names, _ = parse(open(sys.argv[1]).read())
    tips = list(names)
    a, _, _ = far(adj, tips[0])
    b, dist, prev = far(adj, a)
    half = dist[b] / 2
    # walk from b back towards a until the midpoint is passed
    n = b
    while dist[prev[n]] > half:
        n = prev[n]
    p = prev[n]                      # midpoint lies on edge p - n
    l, edge_lab = adj[p][n]
    r = max(dict(adj).keys()) + 1
    del adj[p][n], adj[n][p]
    adj[r][p] = adj[p][r] = [half - dist[p], edge_lab]
    adj[r][n] = adj[n][r] = [dist[n] - half, edge_lab]

    def nwk(node, par):
        if node in names:
            return names[node]
        sub = ",".join(f"{nwk(c, node)}:{adj[node][c][0]:.6g}" if c in names
                       else f"{nwk(c, node)}{adj[node][c][1]}:{adj[node][c][0]:.6g}"
                       for c in adj[node] if c != par)
        return f"({sub})"
    print(nwk(r, None) + ";")


if __name__ == "__main__":
    main()
