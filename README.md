# Formatos de Arquivos em Bioinformática: da referência à variante anotada

Material de aula sobre os principais formatos de arquivo usados em genômica: **FASTA, FASTQ, SAM/BAM/CRAM, VCF/BCF, BED, GFF3 e GTF**. Os comandos são rodados **à mão, passo a passo**, no terminal (bash).

> **Os dados são simulados.** Não vêm de nenhum paciente. As reads foram geradas a partir de um trecho real do genoma humano (hg38), com 114 variantes inseridas de propósito. Assim existe um **gabarito** para conferir se as ferramentas acertaram.

---

## Sumário
1. [Visão geral do fluxo](#1-visão-geral-do-fluxo)
2. [Conceitos e formatos](#2-conceitos-e-formatos)
3. [Os arquivos deste repositório](#3-os-arquivos-deste-repositório)
4. [Ferramentas de cada etapa](#4-ferramentas-de-cada-etapa)
5. [Instalação das ferramentas](#5-instalação-das-ferramentas)
6. [Prática passo a passo](#6-prática-passo-a-passo)
7. [Pontos para discussão](#7-pontos-para-discussão)
8. [Visualizando no IGV](#8-visualizando-no-igv)
9. [Apêndice: como os dados foram gerados](#9-apêndice-como-os-dados-foram-gerados)

---

## 1. Visão geral do fluxo

```
   GENOMA DE REFERÊNCIA          SEQUENCIADOR
        (FASTA)                    (FASTQ)
           │                          │
           │  índices                 │  controle de qualidade
           │  (.fai .dict BWA)        │
           └──────────┬───────────────┘
                      ▼
                 ALINHAMENTO  ──────────►  SAM → BAM → CRAM
                   (BWA)                  (SAMtools)
                      │
                      ▼
            CHAMADA DE VARIANTES  ─────►  VCF / BCF
            (BCFtools, GATK)              (BCFtools)
                      │
                      ▼
                  ANOTAÇÃO  ◄──────────  GFF3 / GTF / BED
               (VEP, ANNOVAR)            (BEDTools)
```

**Região usada:** hg38 **chr6:29.900.001-30.000.000** (100 kb), que contém o gene **HLA-A**. HLA-A é um gene do MHC classe I, importante em transplantes e farmacogenética. No arquivo, a região foi renomeada para `chr6_HLA`, com coordenadas locais: a posição **1** do arquivo corresponde à posição **29.900.001** do chr6.

**Genes na região:**

| Gene | Tipo | Posição local (1-based) | Fita |
|---|---|---|---|
| HCG4B | lncRNA | 24.592-27.215 | − |
| HLA-A | codificante (8 éxons) | 42.532-45.870 | + |
| HCG9 | lncRNA | 75.115-78.403 | + |

---

## 2. Conceitos e formatos

### 2.1 Genoma de referência
O **genoma de referência** é uma sequência montada que representa o genoma de uma espécie. Ela serve de "mapa" para posicionar as reads sequenciadas.

- **Não é o genoma de uma pessoa.** O GRCh38 é um mosaico de DNA de vários doadores (a maior parte de um único indivíduo, o RP11). Ele representa um haplótipo "consenso" e contém alelos raros e lacunas (`N`).
- **Variante = diferença em relação à referência.** Quando um VCF diz `A>G` na posição X, significa que a amostra tem G onde *a referência* tem A. Não quer dizer que G seja "mutante" ou patogênico.
- **Versões (assemblies):** as coordenadas mudam entre versões. **Sempre** informe qual foi usada.

| Nome GRC | Nome UCSC | Ano | Observação |
|---|---|---|---|
| GRCh37 | hg19 | 2009 | Ainda comum em laudos e bancos antigos |
| GRCh38 | hg38 | 2013 | Padrão atual (usado aqui) |
| T2T-CHM13 | hs1 | 2022 | Primeira montagem completa, sem lacunas |

- **Onde baixar:** [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html), [Ensembl](https://www.ensembl.org/info/data/ftp/index.html), [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) e GATK Resource Bundle.
- **Atenção à nomenclatura:** UCSC usa `chr6`, Ensembl usa `6`. Misturar os dois é uma das causas mais comuns de erro ("contig not found").

### 2.2 FASTA (`.fa`, `.fasta`, `.fna`)
Guarda **sequências** (DNA, RNA ou proteína) **sem informação de qualidade**.

```
>chr6_HLA hg38 chr6:29900001-30000000 (regiao HLA-A)      ← cabeçalho: começa com ">"
actattcctccctctgtccctttatatcagttgtttccaggtttgggagtgacattaggt  ← sequência, quebrada em linhas
ttgtctgctgggctggcctagactgcaggcagcaatagtattctagcatgtcttccctca
```
- O primeiro "nome" depois do `>` (até o primeiro espaço) é o **ID** usado por todas as ferramentas (`chr6_HLA`).
- Um FASTA pode ter várias sequências (um genoma tem uma por cromossomo).
- **Minúsculas = soft-masking:** regiões repetitivas (Alu, LINE etc.) marcadas pelo RepeatMasker. A base é a mesma; só fica sinalizado que ali é repetição.
- **`N`** = base desconhecida (lacuna).
- **Índices auxiliares:**
  - `.fai` (`samtools faidx`): permite pular direto para qualquer posição sem ler o arquivo inteiro.
  - `.dict` (`samtools dict` ou `gatk CreateSequenceDictionary`): lista de contigs com tamanho e MD5. É exigido pelo GATK.
  - `.amb .ann .bwt .pac .sa` (`bwa index`): índice FM que o BWA usa para alinhar rápido.

### 2.3 FASTQ (`.fastq`, `.fq`, geralmente `.gz`)
Saída do sequenciador: **sequência + qualidade de cada base**. Cada read ocupa **4 linhas**:

```
@chr6_HLA_hap1-9990/1                  ← 1. identificador (começa com @); /1 = read 1 do par
TAAGAGCCACCCTTAGCTCCTAGAGTTTTATCTC...  ← 2. sequência
+                                      ← 3. separador (pode repetir o ID)
CCCGGGGGGGGGGGCJCJJJJJ1JJJGJJ(JJJG...  ← 4. qualidade: 1 caractere por base
```

**Qualidade Phred:** `Q = −10 × log10(P_erro)`. No padrão **Phred+33**, cada caractere vale `código ASCII − 33`.

| Caractere | ASCII | Q | Probabilidade de erro | Acurácia |
|---|---|---|---|---|
| `+` | 43 | 10 | 1 em 10 | 90% |
| `5` | 53 | 20 | 1 em 100 | 99% |
| `?` | 63 | 30 | 1 em 1.000 | 99,9% |
| `I` | 73 | 40 | 1 em 10.000 | 99,99% |
| `J` | 74 | 41 | — | — |

- **Paired-end:** dois arquivos, `R1` e `R2`. A read N de um arquivo é o par da read N do outro (mesma ordem). As duas vêm das pontas do mesmo fragmento de DNA (aqui, cerca de 400 pb).

### 2.4 SAM / BAM / CRAM
**SAM** (Sequence Alignment/Map) é o formato texto do **alinhamento**: diz *onde* cada read caiu na referência e *como*.

- **BAM:** a versão binária e compactada do SAM. É o formato usado no dia a dia.
- **CRAM:** guarda só as **diferenças em relação à referência**, por isso é ainda menor. Para ler um CRAM, você **precisa do FASTA** usado no alinhamento.

**Cabeçalho** (linhas com `@`):
```
@HD  VN:1.5  SO:coordinate                    ← versão; ordenação (unsorted/coordinate/queryname)
@SQ  SN:chr6_HLA  LN:100000                   ← cada sequência da referência e seu tamanho
@RG  ID:aula  SM:AMOSTRA_SIM  PL:ILLUMINA     ← read group: amostra, plataforma (o GATK exige)
@PG  ID:bwa  PN:bwa  VN:0.7.19  CL:bwa mem... ← programas que geraram/modificaram o arquivo
```

**Corpo:** 11 colunas obrigatórias + tags opcionais:

| # | Coluna | Exemplo | Significado |
|---|---|---|---|
| 1 | QNAME | `chr6_HLA_hap1-9990` | nome da read |
| 2 | FLAG | `83` | soma de bits que descreve a read (ver abaixo) |
| 3 | RNAME | `chr6_HLA` | cromossomo/contig |
| 4 | POS | `96521` | posição (1-based) da 1ª base alinhada |
| 5 | MAPQ | `60` | confiança do mapeamento (Phred; 60 = única, 0 = ambígua) |
| 6 | CIGAR | `150M` | como alinhou: M=match/mismatch, I=inserção, D=deleção, S=soft clip, N=íntron |
| 7 | RNEXT | `=` | cromossomo do par (`=` = mesmo) |
| 8 | PNEXT | `96210` | posição do par |
| 9 | TLEN | `-461` | tamanho do fragmento (negativo = read à direita) |
| 10 | SEQ | `GCCCTGG...` | sequência (já na orientação da referência) |
| 11 | QUAL | `GGCGCC...` | qualidades (Phred+33) |
| — | tags | `NM:i:1 MD:Z:120G29 AS:i:145 RG:Z:aula` | NM = nº de diferenças; MD = onde estão; AS = score; RG = read group |

**FLAG** é uma soma de bits:

| Bit | Valor | Significado |
|---|---|---|
| 0x1 | 1 | read pareada |
| 0x2 | 2 | par alinhado corretamente (*proper pair*) |
| 0x4 | 4 | read **não** mapeada |
| 0x8 | 8 | par não mapeado |
| 0x10 | 16 | read na fita **reversa** |
| 0x20 | 32 | par na fita reversa |
| 0x40 | 64 | é a read 1 |
| 0x80 | 128 | é a read 2 |
| 0x100 | 256 | alinhamento secundário |
| 0x400 | 1024 | duplicata de PCR/óptica |
| 0x800 | 2048 | alinhamento suplementar (quimérico) |

- `83 = 1 + 2 + 16 + 64`: pareada, par correto, fita reversa, read 1.
- `163 = 1 + 2 + 32 + 128`: pareada, par correto, par na reversa, read 2.
- Para decodificar no terminal: `samtools flags 83`, ou acesse <https://broadinstitute.github.io/picard/explain-flags.html>.

### 2.5 VCF / BCF (`.vcf`, `.vcf.gz`, `.bcf`)
**VCF** (Variant Call Format) lista as **posições onde a amostra difere da referência**.

```
##fileformat=VCFv4.2                                    ← meta-informação (##)
##contig=<ID=chr6_HLA,length=100000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM   POS  ID REF ALT QUAL   FILTER INFO          FORMAT       AMOSTRA_SIM   ← cabeçalho (#)
chr6_HLA 1029 .  T   G   222.39 .      DP=29;MQ=60;… GT:PL:DP:AD  0/1:255,0,255:29:14,15
chr6_HLA 1380 .  T   A   225.42 .      DP=34;MQ=60;… GT:PL:DP:AD  1/1:255,102,0:34:0,34
chr6_HLA 3043 .  at  aTGAGTCt 222.06 . …             …            0/1:…
```

| Coluna | Significado |
|---|---|
| CHROM, POS | onde (POS é **1-based**) |
| ID | identificador (ex.: rsID do dbSNP); `.` = nenhum |
| REF, ALT | alelo da referência e alelo alternativo. Em indels, inclui a base anterior (*padding*) |
| QUAL | confiança (Phred) de que existe variante ali |
| FILTER | `PASS`, `.` (não filtrado) ou o nome do filtro que falhou |
| INFO | anotações **do sítio** (DP total, MQ, AF...) |
| FORMAT | **ordem** dos campos da(s) amostra(s) |
| AMOSTRA… | valores **por amostra** |

**Campos de FORMAT mais comuns:**
- **GT** (genótipo): `0/0` = homozigoto referência; `0/1` = heterozigoto; `1/1` = homozigoto alternativo; `./.` = sem chamada. Com `|` (`0|1`) o genótipo está **faseado**.
- **AD** (*allelic depth*): reads com REF,ALT. `14,15` ≈ 50% (heterozigoto); `0,34` = 100% (homozigoto).
- **DP:** profundidade na amostra.
- **PL / GQ:** verossimilhanças dos genótipos e qualidade do genótipo.

**Compressão e índice:**
- `.vcf.gz` deve ser compactado com **bgzip** (não com `gzip` comum) e indexado com **tabix** (`.tbi`) ou `bcftools index` (`.csi`). Assim dá para buscar uma região sem ler o arquivo todo.
- **BCF** é o VCF binário (está para o VCF assim como o BAM está para o SAM).

### 2.6 BED, GFF3 e GTF: anotação e intervalos

| | BED | GFF3 | GTF (GFF2.5) |
|---|---|---|---|
| Coordenadas | **0-based, semiaberto** `[início, fim)` | **1-based, fechado** `[início, fim]` | **1-based, fechado** |
| Uso típico | regiões (alvos de captura, painéis, picos) | anotação hierárquica (gene → mRNA → éxon → CDS) | anotação de genes (Ensembl/GENCODE, RNA-seq) |
| Colunas | 3 obrigatórias (chrom, start, end) + até 12 | 9 fixas | 9 fixas |
| Atributos (col. 9) | — | `ID=...;Parent=...` | `gene_id "..."; transcript_id "...";` |

**O mesmo gene HLA-A nos três formatos:**
```
BED   chr6_HLA  42531  45870  HLA-A  0  +
GFF3  chr6_HLA  RefSeq  gene        42532  45870  .  +  .  ID=gene-HLA-A;Name=HLA-A;gene_biotype=protein_coding
GTF   chr6_HLA  RefSeq  transcript  42532  45870  .  +  .  gene_id "HLA-A"; transcript_id "NM_002116.8";
```
O BED começa em **42531** e o GFF/GTF em **42532**, mas é **exatamente a mesma região**. Para converter: **BED_start = GFF_start − 1** e **BED_end = GFF_end**.

**BED12** acrescenta os blocos (éxons) de um transcrito numa única linha: thickStart/thickEnd (CDS), número de blocos, tamanhos e inícios relativos.

---

## 3. Os arquivos deste repositório

### `FASTA/`: genoma de referência
| Arquivo | Formato | Função | Gerado por | Usado por |
|---|---|---|---|---|
| `chr6_HLA.fa` | FASTA | Sequência de referência (100 kb do chr6, hg38) | UCSC | todas as etapas |
| `chr6_HLA.fa.fai` | índice FAI | Acesso aleatório à sequência: nome, tamanho, offset, bases/linha, bytes/linha | `samtools faidx` | samtools, bcftools, bedtools, GATK, IGV |
| `chr6_HLA.dict` | dicionário SAM | Contigs com tamanho e MD5 | `samtools dict` / `gatk CreateSequenceDictionary` | GATK, Picard |
| `chr6_HLA.fa.amb/.ann/.bwt/.pac/.sa` | índice BWA | Índice FM (Burrows-Wheeler) para alinhamento | `bwa index` | `bwa mem` |

### `FASTQ/`: reads brutas
| Arquivo | Formato | Função |
|---|---|---|
| `amostra_R1.fastq` / `.fastq.gz` | FASTQ (Phred+33) | Read 1 dos pares: 9.990 reads de 150 pb |
| `amostra_R2.fastq` / `.fastq.gz` | FASTQ (Phred+33) | Read 2 dos pares (mesma ordem do R1) |

As versões `.gz` mostram a compressão (3,2 MB → 0,9 MB). As ferramentas leem as duas formas.

### `SAM/`: alinhamento
| Arquivo | Formato | Função |
|---|---|---|
| `amostra.sam` | SAM | Alinhamento em texto, legível (8,1 MB) |
| `amostra.bam` + `.bai` | BAM ordenado + índice | Mesmo conteúdo, binário (1,5 MB). O `.bai` permite buscar regiões |
| `amostra.cram` + `.crai` | CRAM + índice | Compressão por referência (0,7 MB). Precisa do FASTA para abrir |
| `amostra.flagstat.txt` | relatório | Resumo do `samtools flagstat` (100% mapeado) |
| `amostra.depth.tsv` | TSV | Profundidade por posição (`samtools depth -a`), média ≈ 30× |

### `VCF/`: variantes
| Arquivo | Formato | Função |
|---|---|---|
| `variantes.vcf` | VCF 4.2 | Variantes chamadas pelo **bcftools** (114: 91 SNPs e 23 indels) |
| `variantes.vcf.gz` + `.tbi` | VCF bgzip + tabix | Versão compactada e indexada |
| `variantes.bcf` + `.csi` | BCF + índice | VCF binário |
| `variantes.stats.txt` | relatório | Saída do `bcftools stats` |
| `verdade.vcf` / `.vcf.gz` + `.tbi` | VCF | **Gabarito**: as 114 variantes realmente inseridas na simulação |
| `comparacao.txt` | relatório | Chamadas × gabarito, antes e depois da normalização |

### `Anotacao/`: genes da região
| Arquivo | Formato | Função |
|---|---|---|
| `genes.gff3` | GFF3 | Genes, transcritos, éxons e CDS (RefSeq), com hierarquia ID/Parent |
| `genes.gtf` | GTF | Mesma anotação no estilo Ensembl/GENCODE |
| `genes.bed` | BED6 | Um intervalo por gene (0-based) |
| `transcritos.bed12` | BED12 | Transcritos com os éxons em blocos |

### `scripts/`
Usados **apenas para preparar** este material (ver [Apêndice](#9-apêndice-como-os-dados-foram-gerados)). Não fazem parte da aula.

---

## 4. Ferramentas de cada etapa

| Etapa | Ferramenta | O que faz | Entrada → Saída |
|---|---|---|---|
| Indexar referência | **samtools faidx** | Indexa o FASTA e extrai regiões | FASTA → `.fai` |
| | **GATK** `CreateSequenceDictionary` / `samtools dict` | Dicionário de contigs | FASTA → `.dict` |
| | **BWA** `index` | Índice FM para o alinhador | FASTA → `.amb .ann .bwt .pac .sa` |
| Controle de qualidade | FastQC / fastp (citados) | Qualidade por base, adaptadores, duplicação | FASTQ → relatório HTML |
| Alinhamento | **BWA** `mem` | Alinha reads curtas (70 pb–1 Mb) na referência | FASTQ + FASTA → SAM |
| Manipular alinhamentos | **SAMtools** | Converter, ordenar, indexar, filtrar, estatísticas | SAM ↔ BAM ↔ CRAM |
| Marcar duplicatas | **GATK** `MarkDuplicates` (Picard) | Marca cópias de PCR (FLAG 1024) | BAM → BAM |
| Operações com intervalos | **BEDTools** | Interseção, cobertura, extração de sequência, merge | BED/GFF/VCF/BAM → BED/FASTA |
| Chamada de variantes | **BCFtools** `mpileup` + `call` | Chamador estatístico rápido | BAM + FASTA → VCF |
| | **GATK** `HaplotypeCaller` | Remonta haplótipos localmente; padrão-ouro em germinativas | BAM + FASTA → VCF |
| Manipular VCF | **BCFtools** | Filtrar, normalizar, consultar, intersectar, estatísticas | VCF ↔ BCF |
| Filtrar variantes | **GATK** `VariantFiltration` | Marca variantes que falham critérios (FILTER) | VCF → VCF |
| Anotação funcional | **VEP** (Ensembl) | Consequência (missense, sinônima, intrônica…), HGVS, frequências, predições | VCF → VCF/TSV anotado |
| | **ANNOVAR** | Anotação por gene, região e filtro (gnomAD, ClinVar, dbNSFP…) | VCF → tabela (CSV/TXT) |
| Visualização | IGV | Navegador genômico | FASTA + BAM + VCF + GFF |

---

## 5. Instalação das ferramentas

> Nada neste repositório instala coisas automaticamente. Instale só o que for usar.

A forma mais simples (Linux, macOS ou WSL no Windows) é pelo **conda/mamba** com o canal **bioconda**:

```bash
# cria um ambiente separado chamado "formatos"
mamba create -n formatos -c conda-forge -c bioconda samtools bcftools htslib bwa bedtools gatk4
mamba activate formatos

# confere as versões
samtools --version | head -1
bcftools --version | head -1
bwa 2>&1 | grep Version
bedtools --version
gatk --version
```

- **VEP:** `mamba install -c bioconda ensembl-vep`. Para o genoma humano completo é preciso um *cache* (cerca de 20 GB). **Nesta aula não precisa:** usamos o nosso GFF3 (seção 6.6). Também dá para usar a versão web: <https://www.ensembl.org/Tools/VEP>.
- **ANNOVAR:** gratuito para uso acadêmico, mas exige **cadastro** em <https://www.openbioinformatics.org/annovar/annovar_download_form.php>. O link de download chega por e-mail. É um conjunto de scripts Perl (não precisa compilar).
- **IGV:** <https://igv.org/doc/desktop/> (ou a versão web: <https://igv.org/app/>).

---

## 6. Prática passo a passo

> Rode tudo **a partir da pasta raiz do repositório**. Os resultados da prática vão para `pratica/`, assim os arquivos originais ficam intactos para comparação.

```bash
cd Arquivos-de-bioinformatica
mkdir -p pratica
```

### 6.0 Explorando os arquivos com bash puro
Antes de qualquer ferramenta, todo formato aqui é **texto** (exceto BAM/CRAM/BCF) e pode ser inspecionado com comandos básicos.

```bash
# --- FASTA ---
head -n 3 FASTA/chr6_HLA.fa                         # cabeçalho + primeiras linhas
grep -c '^>' FASTA/chr6_HLA.fa                      # quantas sequências? → 1
grep -v '^>' FASTA/chr6_HLA.fa | tr -d '\n' | wc -c # tamanho em bases → 100000

# composição de bases (maiúsculas × minúsculas = soft-masking)
grep -v '^>' FASTA/chr6_HLA.fa | tr -d '\n' | fold -w1 | sort | uniq -c
```
**O que observar:** cerca de 59% das bases estão em minúsculas (repetições mascaradas). A região HLA é rica em elementos repetitivos.

```bash
# --- FASTQ ---
head -n 4 FASTQ/amostra_R1.fastq                      # uma read = 4 linhas
echo $(( $(wc -l < FASTQ/amostra_R1.fastq) / 4 ))    # nº de reads → 9990
zcat FASTQ/amostra_R1.fastq.gz | head -n 8           # ler o .gz sem descompactar
awk 'NR%4==2 {print length($0)}' FASTQ/amostra_R1.fastq | sort | uniq -c   # tamanho das reads → 150

# a read 1 do R1 e a read 1 do R2 são um par (mesmo nome)
head -n 1 FASTQ/amostra_R1.fastq FASTQ/amostra_R2.fastq
```

```bash
# --- SAM / VCF / BED / GFF ---
grep '^@' SAM/amostra.sam                 # cabeçalho do SAM
grep -v '^@' SAM/amostra.sam | head -n 2  # primeiras reads alinhadas
grep -v '^##' VCF/variantes.vcf | head    # VCF sem a meta-informação
column -t Anotacao/genes.bed              # BED alinhado em colunas
awk '$3=="exon"' Anotacao/genes.gff3 | head   # só os éxons do GFF3
```

### 6.1 Preparando a referência: samtools faidx, GATK e BWA

```bash
# Índice FAI
samtools faidx FASTA/chr6_HLA.fa
cat FASTA/chr6_HLA.fa.fai
# chr6_HLA   100000   53   60   61
#   nome     tamanho  offset do 1º byte da sequência   bases/linha   bytes/linha (com \n)

# Extrair uma região (coordenadas 1-based): início do HLA-A
samtools faidx FASTA/chr6_HLA.fa chr6_HLA:42532-42600
# >chr6_HLA:42532-42600
# AGATTCTCCCCAGACGCCGAGGATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACT ...

# Extrair na fita reversa (reverso complementar)
samtools faidx -i FASTA/chr6_HLA.fa chr6_HLA:42532-42600

# Dicionário (escolha UM dos dois)
samtools dict FASTA/chr6_HLA.fa -o FASTA/chr6_HLA.dict
gatk CreateSequenceDictionary -R FASTA/chr6_HLA.fa -O FASTA/chr6_HLA.dict
cat FASTA/chr6_HLA.dict

# Índice do BWA (para o genoma humano inteiro leva cerca de 1 h e uns 5 GB de RAM; aqui, menos de 1 s)
bwa index FASTA/chr6_HLA.fa
ls FASTA/
```

### 6.2 Alinhamento: BWA + SAMtools

```bash
# Alinhar os pares. -R adiciona o read group (obrigatório para o GATK)
bwa mem -t 2 \
  -R '@RG\tID:aula\tSM:AMOSTRA_SIM\tPL:ILLUMINA' \
  FASTA/chr6_HLA.fa FASTQ/amostra_R1.fastq FASTQ/amostra_R2.fastq \
  > pratica/amostra.sam 2> pratica/bwa.log

samtools view -H pratica/amostra.sam              # só o cabeçalho
samtools view pratica/amostra.sam | head -n 2      # só as reads (sem cabeçalho)

# Entendendo a FLAG
samtools flags 83     # 0x53  83  PAIRED,PROPER_PAIR,REVERSE,READ1
samtools flags 163    # 0xa3 163  PAIRED,PROPER_PAIR,MREVERSE,READ2
```

```bash
# SAM → BAM (binário), ordenado por coordenada, e índice
samtools view -b -o pratica/amostra.unsorted.bam pratica/amostra.sam   # só converte
samtools sort -o pratica/amostra.bam pratica/amostra.sam               # converte + ordena
samtools index pratica/amostra.bam                                      # cria .bai
ls -lh pratica/
# amostra.sam 8.1M  |  amostra.unsorted.bam 2.3M  |  amostra.bam 1.5M (ordenado comprime melhor)
```

```bash
# Estatísticas
samtools flagstat pratica/amostra.bam    # 19980 reads, 100% mapeadas, 100% properly paired
samtools idxstats pratica/amostra.bam    # reads por contig (usa o índice, é instantâneo)
samtools coverage pratica/amostra.bam    # cobertura 99,9% da região, profundidade média ≈ 30×
samtools depth -a -r chr6_HLA:42532-42540 pratica/amostra.bam   # profundidade base a base

# Contando e filtrando com samtools view
samtools view -c pratica/amostra.bam             # total de reads      → 19980
samtools view -c -f 16 pratica/amostra.bam       # COM o bit 16 (reversa) → 9990
samtools view -c -F 4 pratica/amostra.bam        # SEM o bit 4 (= mapeadas)
samtools view -c -q 30 pratica/amostra.bam       # MAPQ ≥ 30
samtools view -c pratica/amostra.bam chr6_HLA:42532-45870   # reads no HLA-A → 736 (precisa do .bai)

# Recortar só o HLA-A num BAM novo e voltar para FASTQ
samtools view -b -o pratica/hla_a.bam pratica/amostra.bam chr6_HLA:42532-45870
samtools fastq -1 pratica/hla_R1.fq -2 pratica/hla_R2.fq pratica/hla_a.bam
```

```bash
# BAM → CRAM (precisa da referência para escrever E para ler)
samtools view -C -T FASTA/chr6_HLA.fa -o pratica/amostra.cram pratica/amostra.bam
samtools index pratica/amostra.cram
ls -lh pratica/amostra.sam pratica/amostra.bam pratica/amostra.cram
# 8.1M  →  1.5M  →  727K

# Visualizador de alinhamento no próprio terminal (setas para navegar, "q" para sair)
samtools tview -p chr6_HLA:1020 pratica/amostra.bam FASTA/chr6_HLA.fa
```
**No `tview`:** `.` e `,` são bases iguais à referência (fita + e −), e letras são diferenças. Na posição 1029 aparece metade das reads com `G`: é uma variante heterozigota.

### 6.3 Operações com intervalos: BEDTools

```bash
# Quais variantes caem dentro de genes? (-wa -wb mostra a linha do VCF e a do BED)
bedtools intersect -a VCF/variantes.vcf.gz -b Anotacao/genes.bed -wa -wb | cut -f1,2,4,5,14
# chr6_HLA  42825  G  A  HLA-A
# chr6_HLA  43218  g  gACTAC  HLA-A
# ...

bedtools intersect -a VCF/variantes.vcf.gz -b Anotacao/genes.bed -u | wc -l   # variantes em genes → 11
bedtools intersect -a VCF/variantes.vcf.gz -b Anotacao/genes.bed -v | wc -l   # fora de genes
bedtools intersect -a Anotacao/genes.bed -b VCF/variantes.vcf.gz -c          # nº de variantes por gene
# HCG4B 3 | HLA-A 5 | HCG9 3

# Sequência de cada gene em FASTA (-s respeita a fita: HCG4B sai no reverso complementar)
bedtools getfasta -fi FASTA/chr6_HLA.fa -bed Anotacao/genes.bed -name -s > pratica/genes.fa
head -c 300 pratica/genes.fa

# Cobertura média de cada gene
bedtools coverage -a Anotacao/genes.bed -b pratica/amostra.bam -mean
# HCG4B 32.5 | HLA-A 31.5 | HCG9 30.8

# Cobertura em formato bedGraph (intervalos com a mesma profundidade)
bedtools genomecov -ibam pratica/amostra.bam -bg | head

# Ordenar e juntar intervalos próximos (até 20 kb de distância)
bedtools sort -i Anotacao/genes.bed
bedtools merge -i Anotacao/genes.bed -d 20000

# BEDTools também lê GFF diretamente: unir éxons sobrepostos
awk '$3=="exon"' Anotacao/genes.gff3 | bedtools merge -i - | head
```
**Discussão:** o BEDTools aceita BED (0-based) e GFF/VCF (1-based) **ao mesmo tempo** e converte as coordenadas por dentro. Ele sabe a convenção de cada formato.

### 6.4 Chamada e manipulação de variantes: BCFtools

```bash
# Chamada: mpileup resume as bases por posição; call decide o genótipo
bcftools mpileup -f FASTA/chr6_HLA.fa -a AD,DP pratica/amostra.bam \
  | bcftools call -mv -Oz -o pratica/bcftools.vcf.gz
#   -m = modelo multialélico   -v = só sítios variantes   -Oz = saída VCF compactada (bgzip)
bcftools index -t pratica/bcftools.vcf.gz      # cria .tbi (mesmo que tabix -p vcf)
```

```bash
# Contagens e filtros
bcftools view -H pratica/bcftools.vcf.gz | wc -l                         # total     → 114
bcftools view -H -v snps pratica/bcftools.vcf.gz | wc -l                 # SNPs      → 91
bcftools view -H -v indels pratica/bcftools.vcf.gz | wc -l               # indels    → 23
bcftools view -H -g het pratica/bcftools.vcf.gz | wc -l                  # 0/1       → 67
bcftools view -H -g hom pratica/bcftools.vcf.gz | wc -l                  # 1/1       → 47
bcftools view -H -i 'QUAL>=50 && INFO/DP>=10' pratica/bcftools.vcf.gz | wc -l   # filtro por expressão
bcftools view -H -r chr6_HLA:42532-45870 pratica/bcftools.vcf.gz         # só o HLA-A (usa o índice) → 5

# Transformar o VCF em tabela
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT\t%AD]\n' pratica/bcftools.vcf.gz | head
# chr6_HLA  1029  T  G  222.391  0/1  14,15
# chr6_HLA  1380  T  A  225.417  1/1  0,34

# Estatísticas
bcftools stats pratica/bcftools.vcf.gz | grep '^SN'

# Conversões
bcftools view -Ob -o pratica/bcftools.bcf pratica/bcftools.vcf.gz   # VCF → BCF
bcftools view pratica/bcftools.bcf | head                            # BCF → VCF na tela
```

**Comparando com o gabarito, o problema da representação de indels:**
```bash
bcftools isec -n=2 pratica/bcftools.vcf.gz VCF/verdade.vcf.gz | wc -l    # coincidentes → 106 (?!)

# Normalizar = alinhar indels à esquerda + representação mínima + separar multialélicos
bcftools norm -f FASTA/chr6_HLA.fa -m -any pratica/bcftools.vcf.gz -Oz -o pratica/bcftools.norm.vcf.gz
bcftools index -t pratica/bcftools.norm.vcf.gz
bcftools norm -f FASTA/chr6_HLA.fa -m -any VCF/verdade.vcf.gz -Oz -o pratica/verdade.norm.vcf.gz
bcftools index -t pratica/verdade.norm.vcf.gz

bcftools isec -n=2 pratica/bcftools.norm.vcf.gz pratica/verdade.norm.vcf.gz | wc -l   # → 114 ✔

# Separar em arquivos: só no 1º, só no 2º, nos dois
bcftools isec -p pratica/isec pratica/bcftools.norm.vcf.gz pratica/verdade.norm.vcf.gz
cat pratica/isec/README.txt
```
**Por que 106 → 114?** A mesma inserção pode ser escrita de várias formas. Na posição 3043 o bcftools escreveu `at → aTGAGTCt`, enquanto o gabarito tem `3044 T → TGAGTCT`: é o mesmo evento biológico com outra representação. O `bcftools norm` coloca os dois na forma canônica. **Sempre normalize antes de comparar VCFs ou anotar.**

### 6.5 Chamada de variantes com GATK

```bash
# 1. Marcar duplicatas (reads idênticas vindas da PCR recebem a FLAG 1024)
gatk MarkDuplicates \
  -I pratica/amostra.bam \
  -O pratica/amostra.md.bam \
  -M pratica/md_metrics.txt
samtools index pratica/amostra.md.bam
grep -A2 '^LIBRARY' pratica/md_metrics.txt | cut -f1-9    # 2 pares duplicados (0,02%)

# 2. HaplotypeCaller (precisa de .fai, .dict e do read group no BAM)
gatk HaplotypeCaller \
  -R FASTA/chr6_HLA.fa \
  -I pratica/amostra.md.bam \
  -O pratica/gatk.vcf.gz
bcftools view -H pratica/gatk.vcf.gz | wc -l              # → 114

# 3. Filtro "hard filter": marca a coluna FILTER em vez de apagar a variante
gatk VariantFiltration \
  -R FASTA/chr6_HLA.fa \
  -V pratica/gatk.vcf.gz \
  -O pratica/gatk.filtrado.vcf.gz \
  --filter-name "QD_baixo"   --filter-expression "QD < 2.0" \
  --filter-name "QUAL_baixa" --filter-expression "QUAL < 30.0"
bcftools view -H -f PASS pratica/gatk.filtrado.vcf.gz | wc -l    # quantas passaram

# 4. Comparar GATK × gabarito × BCFtools
bcftools norm -f FASTA/chr6_HLA.fa -m -any pratica/gatk.vcf.gz -Oz -o pratica/gatk.norm.vcf.gz
bcftools index -t pratica/gatk.norm.vcf.gz
bcftools isec -n=2 pratica/gatk.norm.vcf.gz pratica/verdade.norm.vcf.gz   | wc -l   # → 114
bcftools isec -n=2 pratica/gatk.norm.vcf.gz pratica/bcftools.norm.vcf.gz  | wc -l   # → 114
```
**Discussão:** compare as colunas INFO dos dois VCFs (`bcftools view -h` mostra o cabeçalho). O GATK traz QD, FS, SOR, MQRankSum e outras métricas usadas em filtragem. Com dados reais (≈ 30× de cobertura e erros de verdade), os dois chamadores **não** concordam 100%, e em indels o HaplotypeCaller costuma ir melhor.

> Em exomas e genomas inteiros o fluxo GATK Best Practices também inclui **BQSR** (`BaseRecalibrator` + `ApplyBQSR`, com dbSNP e Mills como sítios conhecidos) e, para várias amostras, `-ERC GVCF` + `GenomicsDBImport` + `GenotypeGVCFs`.

### 6.6 Anotação: VEP e ANNOVAR
A chamada diz **onde** está a variante. A anotação diz **o que ela faz**: em qual gene cai, qual o efeito na proteína, a frequência na população e se já foi descrita em doenças.

#### VEP (Ensembl Variant Effect Predictor)
Como nosso genoma é um recorte, usamos o **nosso próprio GFF3** como anotação. O VEP exige o GFF ordenado, compactado com bgzip e indexado:

```bash
# Preparar o GFF3 para o VEP
grep -v '^#' Anotacao/genes.gff3 | sort -k1,1 -k4,4n | bgzip > pratica/genes.gff3.gz
tabix -p gff pratica/genes.gff3.gz

# Anotar (saída VCF: a anotação entra no campo INFO/CSQ)
vep -i pratica/bcftools.norm.vcf.gz \
    --gff pratica/genes.gff3.gz \
    --fasta FASTA/chr6_HLA.fa \
    --vcf --symbol --hgvs --force_overwrite \
    -o pratica/vep.vcf

# Mesma coisa em tabela, com colunas escolhidas
vep -i pratica/bcftools.norm.vcf.gz \
    --gff pratica/genes.gff3.gz \
    --fasta FASTA/chr6_HLA.fa \
    --tab --symbol --hgvs --force_overwrite \
    --fields "Uploaded_variation,Location,Allele,SYMBOL,Consequence,HGVSc,HGVSp" \
    -o pratica/vep.tsv

# Resumo das consequências
grep -v '^#' pratica/vep.tsv | cut -f5 | sort | uniq -c | sort -rn
```
Consequências esperadas: `intergenic_variant`, `intron_variant`, `missense_variant`, `synonymous_variant`, `non_coding_transcript_exon_variant`, `frameshift_variant`... O VEP também gera um relatório HTML (`*_summary.html`).

**Com o genoma humano completo** (dados reais), troque `--gff/--fasta` por `--cache --assembly GRCh38` e adicione, por exemplo, `--af_gnomade --check_existing --sift b --polyphen b`.

#### Coordenadas reais do hg38 (para VEP web e ANNOVAR)
Os bancos públicos usam `chr6` e posições reais. Basta trocar o nome e somar 29.900.000:

```bash
bcftools view pratica/bcftools.norm.vcf.gz \
| awk 'BEGIN{OFS="\t"}
       /^##contig/ {print "##contig=<ID=chr6,length=170805979>"; next}
       /^#/        {print; next}
                   {$1="chr6"; $2=$2+29900000; $4=toupper($4); $5=toupper($5); print}' \
> pratica/variantes_hg38.vcf

grep -v '^#' pratica/variantes_hg38.vcf | head -3 | cut -f1-5
# chr6  29901029  .  T  G
# chr6  29901380  .  T  A
```
O arquivo `pratica/variantes_hg38.vcf` pode ser enviado direto ao **VEP web** (<https://www.ensembl.org/Tools/VEP>, espécie *Human GRCh38*).

#### ANNOVAR
Depois de baixar e descompactar o ANNOVAR (pasta `annovar/`):

```bash
# 1. Baixar os bancos (uma vez só)
perl annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene annovar/humandb/
perl annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20240917 annovar/humandb/
perl annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad41_exome annovar/humandb/

# 2a. Caminho "passo a passo": VCF → formato avinput → anotação por gene
perl annovar/convert2annovar.pl -format vcf4 pratica/variantes_hg38.vcf > pratica/variantes.avinput
head -3 pratica/variantes.avinput        # chr  início  fim  ref  alt  (1-based)
perl annovar/annotate_variation.pl -geneanno -buildver hg38 -dbtype refGene \
     pratica/variantes.avinput annovar/humandb/
# gera pratica/variantes.avinput.variant_function (região: exonic, intronic, intergenic...)
#  e pratica/variantes.avinput.exonic_variant_function (efeito: nonsynonymous SNV, frameshift...)

# 2b. Caminho "tudo de uma vez": várias bases numa tabela só
perl annovar/table_annovar.pl pratica/variantes_hg38.vcf annovar/humandb/ \
     -buildver hg38 -out pratica/annovar \
     -protocol refGene,clinvar_20240917,gnomad41_exome \
     -operation g,f,f \
     -remove -nastring . -vcfinput -polish
# saídas: pratica/annovar.hg38_multianno.txt (tabela) e .hg38_multianno.vcf
```
`-operation`: `g` = anotação por **g**ene, `r` = por **r**egião (ex.: cytoBand), `f` = **f**iltro (bancos de variantes: frequência, ClinVar, predições). Os nomes/versões dos bancos mudam com o tempo; a lista atual está em <https://annovar.openbioinformatics.org/en/latest/user-guide/download/>.

| | VEP | ANNOVAR |
|---|---|---|
| Licença | Livre (Apache 2.0) | Gratuito só para uso acadêmico (com cadastro) |
| Saída | VCF (campo CSQ), TSV, JSON | Tabela TXT/CSV (uma coluna por banco) + VCF |
| Ponto forte | Transcritos Ensembl/RefSeq/MANE, plugins (SpliceAI, CADD, LOFTEE) | Muito prático para juntar dezenas de bancos numa tabela |

---

## 7. Pontos para discussão

1. **Coordenadas 0-based × 1-based:** o HLA-A começa em 42531 no BED e em 42532 no GFF/GTF/VCF/SAM. Erro de ±1 é o bug mais clássico da bioinformática.
2. **Normalização de indels:** sem normalizar, só 106 de 114 variantes batem com o gabarito. Depois do `bcftools norm`, batem as 114 (ver `VCF/comparacao.txt`).
3. **Compressão:** SAM 8,1 MB → BAM 1,5 MB → CRAM 0,7 MB. FASTQ 3,2 MB → `.gz` 0,9 MB. Num genoma humano a 30×, o FASTQ.gz tem cerca de 60–100 GB, o BAM cerca de 60–90 GB e o CRAM cerca de 15–30 GB.
4. **Soft-masking:** as minúsculas do FASTA aparecem até no REF de alguns indels (`at → aTGAGTCt`).
5. **BCFtools × GATK:** os dois acertaram tudo aqui, porque os dados simulados são "limpos". Com dados reais, quais diferenças você esperaria, e por quê?
6. **Nome do cromossomo:** `chr6_HLA` × `chr6` × `6`. O que acontece se o VCF usa `chr6` e o banco de anotação usa `6`?

---

## 8. Visualizando no IGV
1. *Genomes → Load Genome from File* → `FASTA/chr6_HLA.fa`
2. *File → Load from File* → `SAM/amostra.bam`, `VCF/variantes.vcf.gz`, `VCF/verdade.vcf.gz` e `Anotacao/genes.gff3`
3. Na caixa de busca, digite `chr6_HLA:42532-45870` (HLA-A) ou `chr6_HLA:3000-3100` (a inserção do exemplo de normalização).

---

## 9. Apêndice: como os dados foram gerados
Os scripts da pasta `scripts/` serviram **apenas para preparar este material**. Na aula, tudo é feito à mão pelos comandos da seção 6.

1. **Referência:** 100 kb do hg38 (chr6:29.900.001-30.000.000) baixados da API da UCSC. A anotação vem dos genes RefSeq Curated da UCSC, com as coordenadas convertidas para o recorte.
2. **Variantes:** `scripts/simular_variantes.py` inseriu 114 variantes (≈ 80% SNPs, 20% indels de 1 a 6 pb; heterozigotas e homozigotas; seed 42) e gerou dois haplótipos + `verdade.vcf`.
3. **Reads:** `art_illumina` (perfil HiSeq 2500, 2×150 pb, fragmento ≈ 400 pb, ≈ 15× por haplótipo).
4. **Processamento:** `scripts/gerar_dados.sh` (BWA, SAMtools e BCFtools), rodado num container Docker definido em `scripts/Dockerfile`.
