#!/bin/bash
# crispr-geno launcher for macOS.
# Double-click this file in Finder to run the tool in a Terminal window.
#
# First-time users: macOS will warn the file is "from an unidentified developer".
# Right-click the file -> Open -> confirm. After that, double-click works.

set -e

# Activate the dedicated conda env
CONDA_BASE="/opt/homebrew/Caskroom/miniconda/base"
if [ ! -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    # Fallback for Intel Macs with /usr/local miniconda
    CONDA_BASE="$HOME/miniconda3"
fi
if [ ! -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    echo "ERROR: could not find conda. Edit this file and set CONDA_BASE to your miniconda/anaconda path."
    echo "Press Return to close."
    read _
    exit 1
fi

source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate crispr-geno

clear
cat <<'BANNER'

  ╔══════════════════════════════════════════════╗
  ║                 crispr-geno                  ║
  ║  Long-read CRISPR amplicon genotyping tool   ║
  ╚══════════════════════════════════════════════╝

  What it does
    • Reads Nanopore long-read amplicon fastqs
    • Aligns them to your amplicon with minimap2
    • Calls per-plant, per-guide allele genotypes
    • Flags plants that still carry active Cas9

  How to use it
    • Put a samples sheet (samples.txt / .tsv / .csv / .xlsx)
      in the folder with your fastqs
    • Columns: name, fastq (or fastq_r1), amplicon_seq, guide_seq
      (guide_seq is comma-separated)
    • Results go to <your-folder>/CRISPRgeno_results/

  Willian Viana
  Dinneny Lab
  contact: williangviana@outlook.com

BANNER

# Three ways to give the tool a folder:
#   (a) drop the folder onto this .command file in Finder (first run only; $1),
#   (b) paste / drag-and-drop the folder into this Terminal window,
#   (c) just press Enter to open a native macOS folder picker.
# After each run you can process another folder without closing the window.

FIRST_RUN=1

while true; do
    INPUT_DIR=""

    # (a) folder dropped onto the .command file — only honored on the first pass
    if [ "$FIRST_RUN" = "1" ] && [ $# -ge 1 ] && [ -n "$1" ] && [ -d "$1" ]; then
        INPUT_DIR="$1"
    fi
    FIRST_RUN=0

    if [ -z "$INPUT_DIR" ]; then
        echo "Drag your fastq folder into this window, then press Enter."
        echo "  (or just press Enter to open a folder picker, or type 'q' to quit)"
        read -r DROPPED

        # Quit
        if [ "$DROPPED" = "q" ] || [ "$DROPPED" = "Q" ] || [ "$DROPPED" = "quit" ]; then
            break
        fi

        # Trim quotes, un-escape spaces, strip trailing whitespace
        DROPPED="${DROPPED%\'}"; DROPPED="${DROPPED#\'}"
        DROPPED="${DROPPED%\"}"; DROPPED="${DROPPED#\"}"
        DROPPED="${DROPPED//\\ / }"
        DROPPED="${DROPPED%"${DROPPED##*[![:space:]]}"}"

        if [ -n "$DROPPED" ] && [ -d "$DROPPED" ]; then
            INPUT_DIR="$DROPPED"
        elif [ -n "$DROPPED" ]; then
            echo "Not a folder: $DROPPED"
            echo "Falling back to folder picker..."
        fi
    fi

    # (c) folder picker
    if [ -z "$INPUT_DIR" ]; then
        INPUT_DIR=$(osascript -e 'POSIX path of (choose folder with prompt "Select the folder containing your fastq files")' 2>/dev/null || true)
    fi

    if [ -z "$INPUT_DIR" ] || [ ! -d "$INPUT_DIR" ]; then
        echo "No folder selected."
        echo
        continue
    fi
    INPUT_DIR="${INPUT_DIR%/}"

    OUTPUT_DIR="${INPUT_DIR}/CRISPRgeno_results"

    echo
    echo "Input  : $INPUT_DIR"
    echo "Output : $OUTPUT_DIR"
    echo

    crispr-geno --input-dir "$INPUT_DIR" -o "$OUTPUT_DIR" || true

    echo
    echo "Opening results folder..."
    open "$OUTPUT_DIR" 2>/dev/null || true

    echo
    echo "---"
done

echo
echo "Press Return to close this window."
read _
