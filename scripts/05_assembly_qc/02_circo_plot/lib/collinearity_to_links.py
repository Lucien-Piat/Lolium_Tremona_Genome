#!/usr/bin/env python3
"""Convert an MCScanX .collinearity file to a Circos links TSV."""
import re
import sys

HEADER_RE = re.compile(r"^## Alignment\s+\d+:.*?(plus|minus)\s*$")

def load_gene_positions(gff_path):
    pos = {}
    with open(gff_path) as fh:
        for line in fh:
            p = line.rstrip().split("\t")
            if len(p) >= 4:
                pos[p[1]] = (p[0], int(p[2]), int(p[3]))
    return pos

def parse_collinearity(path, gene_pos):
    blocks = []
    current = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith("###"):
                continue
            if line.startswith("##"):
                m = HEADER_RE.match(line)
                current = {"orient": m.group(1) if m else "plus", "pairs": []}
                blocks.append(current)
                continue
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3 and current is not None:
                g1, g2 = parts[1].strip(), parts[2].strip()
                if g1 in gene_pos and g2 in gene_pos:
                    current["pairs"].append((g1, g2))
    return blocks

def block_to_link(block, gene_pos, prefix):
    if not block["pairs"]:
        return None
    c1s, c2s = set(), set()
    s1_min = s2_min = float("inf")
    e1_max = e2_max = 0
    for g1, g2 in block["pairs"]:
        c1, s1, e1 = gene_pos[g1]
        c2, s2, e2 = gene_pos[g2]
        c1s.add(c1); c2s.add(c2)
        s1_min, e1_max = min(s1_min, s1), max(e1_max, e1)
        s2_min, e2_max = min(s2_min, s2), max(e2_max, e2)
    if len(c1s) != 1 or len(c2s) != 1:
        return None
    c1 = re.sub(rf"^{prefix}", "chr", c1s.pop())
    c2 = re.sub(rf"^{prefix}", "chr", c2s.pop())
    strand = "+" if block["orient"] == "plus" else "-"
    return (c1, s1_min, e1_max, c2, s2_min, e2_max, strand, len(block["pairs"]))

def main(name, prefix, out_path):
    gene_pos = load_gene_positions(f"{name}.gff")
    blocks = parse_collinearity(f"{name}.collinearity", gene_pos)
    n_total = n_intra = 0
    with open(out_path, "w") as out:
        out.write("q_chr\tq_start\tq_end\tt_chr\tt_start\tt_end\tstrand\tn_genes\n")
        for b in blocks:
            link = block_to_link(b, gene_pos, prefix)
            if link is None:
                continue
            if link[0] == link[3]:
                n_intra += 1
            out.write("\t".join(str(x) for x in link) + "\n")
            n_total += 1
    print(f"Wrote {n_total} blocks ({n_intra} intra-chromosome) to {out_path}",
          file=sys.stderr)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])