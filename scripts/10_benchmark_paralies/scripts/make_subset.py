target_bp = int(float(snakemake.params.target_mb) * 1e6)

sequences, order = {}, []
name, chunks = None, []
with open(snakemake.input.fa) as fh:
    for line in fh:
        if line.startswith(">"):
            if name is not None:
                sequences[name] = "".join(chunks)
            name = line[1:].split()[0]
            order.append(name)
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        sequences[name] = "".join(chunks)

keep = {}
used = 0
for n in order:
    if used >= target_bp:
        break
    room = target_bp - used
    take = min(len(sequences[n]), room)
    if take < 1_000_000:
        break
    keep[n] = take
    used += take

if not keep:
    raise SystemExit(f"target {target_bp:,} bp is smaller than one chromosome")


def wrap(fh, name, seq, width=60):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


with open(snakemake.output.fa, "w") as out:
    for n in order:
        if n in keep:
            wrap(out, n, sequences[n][:keep[n]])

with open(snakemake.input.gff) as fin, open(snakemake.output.gff, "w") as out:
    for line in fin:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) != 9 or p[0] not in keep:
            continue
        if int(p[4]) > keep[p[0]]:
            continue
        out.write(line)

