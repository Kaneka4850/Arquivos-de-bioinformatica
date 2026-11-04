#!/usr/bin/env bash
# Roda DENTRO do container (diretorio Formatos montado em /data).
set -euo pipefail
cd /data
REF=FASTA/chr6_HLA.fa
# 1) Indices da referencia
samtools faidx $REF
samtools dict $REF -o FASTA/chr6_HLA.dict
bwa index $REF 2> /dev/null
# 2) Variantes conhecidas -> haplotipos + gabarito
mkdir -p tmp
python scripts/simular_variantes.py $REF tmp
cp tmp/verdade.vcf VCF/verdade.vcf
# 3) Reads Illumina paired-end 2x150, ~15x por haplotipo (~30x total)
for h in 1 2; do
  art_illumina -ss HS25 -i tmp/hap$h.fa -p -l 150 -f 15 -m 400 -s 50 -rs $((41+h)) -na -q -o tmp/h${h}_ > /dev/null
done
cat tmp/h1_1.fq tmp/h2_1.fq > FASTQ/amostra_R1.fastq
cat tmp/h1_2.fq tmp/h2_2.fq > FASTQ/amostra_R2.fastq
gzip -kf FASTQ/amostra_R1.fastq FASTQ/amostra_R2.fastq
# 4) Alinhamento: SAM -> BAM ordenado + indice -> CRAM
bwa mem -t 2 -R '@RG\tID:sim\tSM:AMOSTRA_SIM\tPL:ILLUMINA' $REF FASTQ/amostra_R1.fastq FASTQ/amostra_R2.fastq 2> /dev/null > SAM/amostra.sam
samtools sort -o SAM/amostra.bam SAM/amostra.sam
samtools index SAM/amostra.bam
samtools view -C -T $REF -o SAM/amostra.cram SAM/amostra.bam
samtools index SAM/amostra.cram
samtools flagstat SAM/amostra.bam > SAM/amostra.flagstat.txt
samtools depth -a SAM/amostra.bam > SAM/amostra.depth.tsv
# 5) Chamada de variantes: VCF, VCF.gz + tabix, BCF
bcftools mpileup -f $REF -a AD,DP SAM/amostra.bam 2> /dev/null | bcftools call -mv -Ov -o VCF/variantes.vcf
bgzip -kf VCF/variantes.vcf && tabix -fp vcf VCF/variantes.vcf.gz
bcftools view -Ob -o VCF/variantes.bcf VCF/variantes.vcf.gz && bcftools index -f VCF/variantes.bcf
bcftools stats VCF/variantes.vcf.gz > VCF/variantes.stats.txt
# 6) Comparacao com o gabarito
bgzip -kf VCF/verdade.vcf && tabix -fp vcf VCF/verdade.vcf.gz
bcftools norm -f $REF -m -any VCF/variantes.vcf.gz -Oz -o tmp/var.norm.vcf.gz 2> /dev/null && tabix -p vcf tmp/var.norm.vcf.gz
bcftools norm -f $REF -m -any VCF/verdade.vcf.gz -Oz -o tmp/ver.norm.vcf.gz 2> /dev/null && tabix -p vcf tmp/ver.norm.vcf.gz
n=$(grep -vc '^#' VCF/verdade.vcf); c=$(grep -vc '^#' VCF/variantes.vcf)
bruto=$(bcftools isec -n=2 VCF/variantes.vcf.gz VCF/verdade.vcf.gz | wc -l)
norm=$(bcftools isec -n=2 tmp/var.norm.vcf.gz tmp/ver.norm.vcf.gz | wc -l)
gt=$(bcftools isec -n=2 -w1 tmp/var.norm.vcf.gz tmp/ver.norm.vcf.gz -Oz -o tmp/a.vcf.gz && tabix -p vcf tmp/a.vcf.gz &&      paste <(bcftools query -f '%POS %REF %ALT [%GT]
' tmp/a.vcf.gz)            <(bcftools query -R tmp/a.vcf.gz -f '%POS %REF %ALT [%GT]
' tmp/ver.norm.vcf.gz) | awk '$4==$8' | wc -l)
{ echo "Gabarito (verdade.vcf):            $n variantes"
  echo "Chamadas pelo bcftools:            $c variantes"
  echo "Coincidentes SEM normalizar:       $bruto  (indels escritos de forma diferente)"
  echo "Coincidentes APOS bcftools norm:   $norm"
  echo "Genotipo (0/1 x 1/1) correto:      $gt"; } | tee VCF/comparacao.txt
rm -rf tmp
