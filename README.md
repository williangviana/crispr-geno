# crispr-geno

Genotype long-read CRISPR amplicon sequencing (Nanopore / PacBio). For each sample, phases reads into the underlying chromosomes, calls per-chromosome edits and large structural deletions, and writes a single homozygous / biallelic / heterozygous / WT zygosity call per sample alongside the per-guide detail in a color-coded XLSX.

Built for amplicons too long for CRISPResso2 (> ~2 kb). Alignment via minimap2, so reads from ~200 bp up to the full amplicon work.

## Install

macOS one-liner — sets up Homebrew, Python 3.12, minimap2, samtools, and the `crispr-geno` command:

```bash
curl -fsSL https://raw.githubusercontent.com/williangviana/crispr-geno/stable/install/install.sh | bash
```

Open a **new** Terminal after install, then `crispr-geno --help`.

<details>
<summary>Conda alternative</summary>

```bash
git clone https://github.com/williangviana/crispr-geno.git
cd crispr-geno
# Apple Silicon: prefix with CONDA_SUBDIR=osx-64 (runs under Rosetta)
CONDA_SUBDIR=osx-64 conda env create -f environment.yml
conda activate crispr-geno
```
</details>

## Usage

Drop a tab-separated `samples.txt` into the folder with your fastqs:

```
name	fastq	amplicon_seq	guide_seq
plant-1	plant-1.fastq	CTGCATAGC...TGTGG	AAGCTATTCTCCCTCGGGGT,TCACACGAGCGTAACAAGGA
plant-2	plant-2.fastq	CTGCATAGC...TGTGG	AAGCTATTCTCCCTCGGGGT,TCACACGAGCGTAACAAGGA
```

- The fastq column can be `fastq` or `fastq_r1` (CRISPResso batch sheets work unchanged).
- `guide_seq` is comma-separated.
- Amplicon + guides are taken from the first row — all samples share them.

Then double-click `crispr-geno.command` and pick the folder, or:

```bash
crispr-geno --input-dir /path/to/fastq_folder -o results/
```

### Other modes

- **No samples.txt** — the CLI prompts for amplicon + guides interactively.
- **Non-interactive flags:**
  ```bash
  crispr-geno \
      --amplicon ACGT...TTT \
      --guide gRNA1:AAGCTATTCTCCCTCGGGGT \
      --guide gRNA2:TCACACGAGCGTAACAAGGA \
      --input-dir "T2 sequencing/" -o T2_results/
  ```
  (`--amplicon` also accepts a FASTA path.)
- **Scattered fastqs** — TSV mapping `sample_id<TAB>fastq_path`:
  ```bash
  crispr-geno --samples-tsv samples.tsv --input-dir ignored -o results/
  ```

### Tuning

| Flag | Default | Purpose |
|---|---:|---|
| `--window` | 5 | bp each side of the cut where indels count |
| `--min-spanning` | 10 | min reads spanning the cut to make a call |
| `--homo-threshold` | 70 | % required for a homozygous call |
| `--het-threshold` | 20 | % required for a het call |
| `--min-read-len` | 200 | drop reads shorter than this |
| `--threads` | 4 | threads for minimap2 / samtools |
| `--phasing-noise-frac` | 0.05 | merge haplotype tuples below this fraction into the nearest dominant tuple |
| `--phasing-noise-min-reads` | 2 | minimum reads for a tuple to escape noise merging |
| `--mosaic-threshold` | 0.10 | flag sample as mosaic when ≥3 haplotype tuples each cross this fraction |
| `--legacy-output` | (off) | make the legacy per-guide table the primary xlsx sheet (compat for old downstream consumers) |

## Outputs

In the `--output` directory:

| File | Contents |
|---|---|
| `results.xlsx` | Primary "Sequences" sheet: per-sample zygosity call (`homozygous` / `biallelic` / `heterozygous` / `WT` / `mosaic` / `lowN`), per-chromosome edit summaries (`Allele 1`, `Allele 2`), and per-guide `A1 / A2` cells. The legacy per-guide table is preserved as a hidden "Per-guide raw" sheet. |
| `editing_rates.csv` | Per-sample × per-guide indel counts, the legacy per-guide call, supporting read fraction, and `call_confidence` |
| `sv_events.csv` | One row per detected large deletion (≥50 bp) spanning two cut sites, with supporting / spanning read counts and the SV call. |
| `haplotypes.csv` | One row per inferred haplotype tuple per sample (post noise-merging), with read support and percent. The audit trail behind the `Allele 1` / `Allele 2` columns. |
| `reference.fasta` | Amplicon in FASTA |
| `bam/<sample>.bam` | Alignments — cached for reruns |

### Call confidence

| Label | Meaning |
|---|---|
| `high` | ≥100 spanning reads, dominant allele ≥90% |
| `med` | between the low and high thresholds |
| `low` | <50 spanning reads, **or** call barely crossed threshold |
| `ambiguous` | competing allele at ≥15% — may be masking biology (red border in XLSX) |
| `lowN` | fewer than `--min-spanning` reads — no call |

`ambiguous` is the cell to look at: mosaic, chimeric, biallelic-with-minor-WT, or messy data.

### Notation

**Zygosity call (primary results.xlsx column):**

This is a DNA-level statement only — it says nothing about whether the protein is functionally disrupted. Whether an in-frame indel is hypomorphic, neutral, or null is for a downstream protein-level assay to decide.

| Code | Meaning |
|---|---|
| `homozygous` | both chromosomes carry the *same* edit — breeds true |
| `biallelic` | both chromosomes are edited but carry *different* lesions (segregates in progeny) |
| `heterozygous` | one chromosome edited, the other wildtype |
| `WT` | both chromosomes wildtype — no edits detected |
| `mosaic` | ≥3 distinct haplotypes — Cas9 still active or chimeric tissue |
| `lowN` | fewer than `--min-spanning` reads phased |

**Per-chromosome `Allele 1` / `Allele 2`:**
- `WT` — no edits on this chromosome
- `SV:-Nbp` — large deletion of N bp spanning two guides
- `gRNAi:+SEQ` / `gRNAi:-SEQ` — insertion / deletion of `SEQ` at the i-th guide cut
- comma-joined when one chromosome carries edits at multiple guides

**Per-guide `A1 / A2` cells** show the slot label on each chromosome, separated by `/` (e.g. `SV / WT`, `SV / +T`, `WT / WT`).

**Legacy per-guide call** (hidden "Per-guide raw" sheet, also in `editing_rates.csv`):

| Code | Meaning |
|---|---|
| `WT` | no edit above noise |
| `+XYZ` / `-XYZ` | homozygous insertion / deletion of `XYZ` at the cut |
| `+XYZ/WT` / `-XYZ/WT` | heterozygous — one chromosome edited, one wildtype |
| `SV:-Nbp` / `SV:-Nbp/WT` | structural variant column for a guide pair |
| `lowN` | too few reads for a call |

## Troubleshooting

- **`guide not found in amplicon`** — check the guide sequence; both strands are searched.
- **`minimap2 not found on PATH`** — reinstall via the one-liner, or activate the conda env.
- **Very low spanning counts** — reads may be too short. Try lowering `--min-read-len`.
- **A gRNA shows no edits anywhere** — it may not have cut in your parental line, or the edit falls outside the amplicon.
- **Genotype says `del` but sequence shows `WT` or `/+`** — many distinct deletions summing to a high aggregate but no single dominant allele. Signature of a mosaic (Cas9 still active) or a chimera-rich library.

## Out of scope

No translocation or inter-chromosomal calls, no base-editing-specific calls, no multi-amplicon multiplexing in one run. Only diploid samples are modelled; polyploid haplotype recovery is not supported.
