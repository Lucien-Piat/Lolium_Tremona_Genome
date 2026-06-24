#!/usr/bin/env python3
# D-GENIES style dotplot

import argparse # type: ignore
import matplotlib # type: ignore
matplotlib.use("Agg")
import matplotlib.pyplot as plt # type: ignore


def read_fai_order(path, min_len):
    order, length = [], {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            name, n = f[0], int(f[1])
            if n >= min_len:
                order.append(name)
                length[name] = n
    return order, length


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf", required=True)
    ap.add_argument("--qfai", required=True)
    ap.add_argument("--rfai", required=True)
    ap.add_argument("--min-chrom-len", type=int, default=50000000)
    ap.add_argument("--min-block", type=int, default=300000)
    ap.add_argument("--title", default="Tremona vs reference")
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    qorder, qlen = read_fai_order(args.qfai, args.min_chrom_len)
    rorder, rlen = read_fai_order(args.rfai, args.min_chrom_len)

    qpair = dict(zip(qorder, rorder))   # query chrom -> fai-paired reference chrom
    qkeep = set(qorder)

    seg = {}      # q -> list of (qs, qe, ts, te, strand)
    with open(args.paf) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 9:
                continue
            q, qs, qe, st, t = f[0], int(f[2]), int(f[3]), f[4], f[5]
            ts, te = int(f[7]), int(f[8])
            if q not in qkeep or (qe - qs) < args.min_block:
                continue
            if t != qpair.get(q):
                continue
            seg.setdefault(q, []).append((qs, qe, ts, te, st))

    ncol = 3
    nrow = (len(qorder) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 5 * nrow))
    axes = axes.ravel() if hasattr(axes, "ravel") else [axes]

    summary = []
    for ax, q in zip(axes, qorder):
        r = qpair[q]
        blocks = seg.get(q, [])
        plus = sum(qe - qs for qs, qe, ts, te, st in blocks if st == "+")
        minus = sum(qe - qs for qs, qe, ts, te, st in blocks if st == "-")
        orient = "forward" if plus >= minus else "INVERTED"
        summary.append((q, r, orient, plus, minus))

        if not blocks:
            ax.set_title("%s vs %s : no blocks" % (q, r), fontsize=9)
            ax.tick_params(labelsize=6)
            continue
        for qs, qe, ts, te, st in blocks:
            if st == "+":
                ax.plot([ts, te], [qs, qe], color="0.15", lw=1.4)
            else:
                ax.plot([ts, te], [qe, qs], color="0.15", lw=1.4)
        ax.set_xlim(0, rlen[r]); ax.set_ylim(0, qlen[q])
        ax.set_title("%s vs %s" % (q, r), fontsize=9)
        ax.set_xlabel("%s pos (bp)" % r, fontsize=8)
        ax.set_ylabel("%s pos (bp)" % q, fontsize=8)
        ax.tick_params(labelsize=6)

    for ax in axes[len(qorder):]:
        ax.axis("off")

    fig.suptitle("%s, paired by fai order" % args.title, fontsize=13)
    fig.tight_layout()
    fig.savefig(args.out_prefix + ".perchrom.png", dpi=200)
    fig.savefig(args.out_prefix + ".perchrom.pdf")
    plt.close(fig)

    with open(args.out_prefix + ".tsv", "w") as out:
        out.write("query\tref\torientation\tplus_bp\tminus_bp\n")
        for q, r, o, p, m in summary:
            out.write("%s\t%s\t%s\t%d\t%d\n" % (q, r, o, p, m))


if __name__ == "__main__":
    main()