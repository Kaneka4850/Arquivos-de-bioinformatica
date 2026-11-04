"""Insere variantes conhecidas (SNPs e indels, het/hom) em uma referencia
e gera dois haplotipos (hap1.fa, hap2.fa) + o gabarito verdade.vcf."""
import random, sys
ref_fa, outdir = sys.argv[1], sys.argv[2]
random.seed(42)
name, seq = None, []
for line in open(ref_fa):
    if line.startswith(">"): name = line[1:].split()[0]
    else: seq.append(line.strip())
ref = "".join(seq); L = len(ref)
# posicoes espacadas (>=300 pb) longe das bordas
pos, p = [], 500
while True:
    p += random.randint(300, 1500)
    if p > L - 500: break
    pos.append(p)
vars_ = []
for p in pos:  # p = indice 0-based
    r = ref[p].upper()
    if r not in "ACGT": continue
    kind = random.random()
    gt = random.choice(["0/1", "0/1", "1/1"])
    if kind < 0.8:                                   # SNP
        alt = random.choice([b for b in "ACGT" if b != r]); R = r
    elif kind < 0.9:                                 # insercao
        R = r; alt = r + "".join(random.choice("ACGT") for _ in range(random.randint(1, 6)))
    else:                                            # delecao
        n = random.randint(1, 6); R = ref[p:p+n+1].upper(); alt = r
        if not set(R) <= set("ACGT"): continue
    vars_.append((p, R, alt, gt))
def build(hap):
    out, last = [], 0
    for p, R, alt, gt in vars_:
        if gt == "1/1" or (gt == "0/1" and hap == (p % 2)):
            out.append(ref[last:p]); out.append(alt); last = p + len(R)
    out.append(ref[last:]); return "".join(out)
for h in (0, 1):
    s = build(h)
    with open(f"{outdir}/hap{h+1}.fa", "w", newline="\n") as f:
        f.write(f">{name}_hap{h+1}\n" + "\n".join(s[i:i+60] for i in range(0, len(s), 60)) + "\n")
with open(f"{outdir}/verdade.vcf", "w", newline="\n") as f:
    f.write("##fileformat=VCFv4.2\n##source=simular_variantes.py (gabarito)\n")
    f.write(f"##contig=<ID={name},length={L}>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAMOSTRA_SIM\n")
    for p, R, alt, gt in vars_:
        f.write(f"{name}\t{p+1}\tsim{p+1}\t{R}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")
print(f"{len(vars_)} variantes simuladas")
