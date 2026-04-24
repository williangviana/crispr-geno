# crispr-geno

Genotype long-read CRISPR amplicon sequencing (Nanopore / PacBio).

For each plant or sample, finds the dominant indel at each guide's cut site
and outputs:

- per-site indel rates (TSV)
- genotype calls (TSV + color-coded XLSX)
- dominant allele sequences (TSV + color-coded XLSX)

Designed for long-read amplicon sequencing where the amplicon is longer than
tools like CRISPResso2 can handle (> ~2 kb). Uses minimap2 for alignment, so
reads can be any length from ~200 bp up to full amplicon + overhangs.

## Install

Requires `minimap2` and `samtools` on PATH. The repo ships an `environment.yml`
that sets up everything in a dedicated conda env:

```bash
git clone <this-repo-or-copy-the-folder>
cd crispr-geno

# On Apple Silicon, prefix with CONDA_SUBDIR=osx-64 (some bioconda packages
# don't have native arm64 builds yet — runs under Rosetta):
CONDA_SUBDIR=osx-64 conda env create -f environment.yml

conda activate crispr-geno
```

The env name is `crispr-geno`. `pip install -e .` runs automatically as part
of env creation, so the `crispr-geno` command is on PATH as soon as you
activate.

Verify:

```bash
crispr-geno --version
```

## Usage

### Easiest: drop a `samples.txt` in your fastq folder

Create a tab-separated `samples.txt` in the folder that holds your fastqs, with
these four columns:

```
name	fastq_r1	amplicon_seq	guide_seq
plant-1	plant-1.fastq	CTGCATAGC...TGTGG	AAGCTATTCTCCCTCGGGGT,TCACACGAGCGTAACAAGGA
plant-2	plant-2.fastq	CTGCATAGC...TGTGG	AAGCTATTCTCCCTCGGGGT,TCACACGAGCGTAACAAGGA
...
```

Guide sequences are comma-separated. The amplicon and guides are taken from the
first row (all samples must share the same amplicon and guides in a single run).
This is the same format CRISPResso's batch mode uses, so if you already have
one of those you can drop it in unchanged.

Then either double-click `crispr-geno.command` and pick the folder, or run:

```bash
crispr-geno --input-dir /path/to/fastq_folder -o results/
```

No prompts — it reads everything from `samples.txt`.

### Or: fully interactive

If there's no `samples.txt` in the folder, the CLI prompts you to paste:

1. The amplicon sequence (one or more lines; blank line to finish).
2. Your guide sequences (`NAME SEQ` one per line; blank line to finish).

### Fully non-interactive

```bash
crispr-geno \
    --amplicon ACGTACGT...TTT \
    --guide gRNA1:AAGCTATTCTCCCTCGGGGT \
    --guide gRNA2:TCACACGAGCGTAACAAGGA \
    --guide gRNA3:GCCAGCTTTCAACCTGTTCC \
    --guide gRNA4:GATGGGAGCAACCCGAGTAC \
    --input-dir  "T2 sequencing/" \
    -o T2_results/
```

You can also pass the amplicon as a FASTA file path.

### Sample list via TSV

If your fastq files are scattered, point at a TSV instead of a directory:

```
sample_id<TAB>fastq_path
1-296-1<TAB>/data/run1/1-296-1.fastq
2-296-2<TAB>/data/run2/2-296-2.fastq
...
```

Then:

```bash
crispr-geno --samples-tsv samples.tsv --input-dir ignored -o results/
```

### Tuning

| Flag | Default | What it does |
|---|---:|---|
| `--window` | 5 | bp each side of the cut site where indels "count" |
| `--min-spanning` | 20 | min reads spanning the cut to make a call |
| `--homo-threshold` | 70 | % required for a homozygous call |
| `--het-threshold` | 20 | % required for a het call |
| `--min-read-len` | 200 | drop reads shorter than this (adapter dimers) |
| `--threads` | 4 | threads for minimap2 / samtools |

## Outputs

All in the `--output` directory:

| File | Contents |
|---|---|
| `editing_rates.csv` | Per-sample × per-guide indel counts, plus the final `call`, the fraction of reads supporting it, and a `call_confidence` label (`high`/`med`/`low`/`ambiguous`/`lowN`) |
| `results.xlsx` | Color-coded dominant allele sequence at each site (`-GATC`, `+T`, etc.). Includes an inline Legend section and a red border around any cell flagged `ambiguous`. |
| `reference.fasta` | The amplicon written to FASTA for re-use |
| `bam/<sample>.bam` | Alignments — cached for reruns |

### Call confidence

Each call in `editing_rates.csv` gets a confidence label:

| Label | Meaning |
|---|---|
| `high` | ≥100 spanning reads and dominant allele ≥90% of them |
| `med` | most cases between the low/high thresholds |
| `low` | <50 spanning reads, OR call barely crossed threshold (70–80% for homo) |
| `ambiguous` | a competing allele is present at ≥15% — the call may be masking biology. These cells also get a red border in the xlsx |
| `lowN` | fewer than `--min-spanning` reads — no call possible |

An `ambiguous` cell is the one to look at: either the plant is more complex than a single-allele call captures (mosaic, chimeric, or biallelic with a minor WT population), or the data is messy.

## Notation in the tables

| Code | Meaning |
|---|---|
| `WT` | no edit above noise |
| `+XYZ` | homozygous insertion of bases `XYZ` at the cut |
| `+XYZ/WT` | heterozygous: one chromosome has +XYZ, the other is wildtype |
| `-XYZ` | homozygous deletion of bases `XYZ` from the reference |
| `-XYZ/WT` | heterozygous: one chromosome has -XYZ, the other is wildtype |
| `lowN` | fewer than `--min-spanning` reads |

## When a genotype call and a sequence call disagree

If `genotype_table` says `del` but `sequence_table` shows `WT` or `/+`, the
sample has *many* different deletion indels summing to a high aggregate rate,
but no single dominant allele. This is the signature of a mosaic plant with
Cas9 still active, or a PCR chimera-rich library.

## Troubleshooting

- **"guide not found in amplicon"**: the guide you provided isn't in your
  amplicon on either strand. Double-check the sequence.
- **"minimap2 not found on PATH"**: activate the conda env.
- **Very low spanning counts** (<20 for most samples): your reads may be too
  short. Try lowering `--min-read-len`, but remember short reads at a long
  amplicon often can't span all cut sites.
- **gRNA shows no edits in any sample**: that guide may simply not have cut
  in your parental line, or the edit may fall outside the amplicon.

## What it doesn't do

- No structural-variant or translocation detection. For large deletions
  spanning multiple cut sites, use CRISPRlungo or similar.
- No base-editing specific calls (only indels around the cut).
- No multiplexing across multiple amplicons in the same run.
