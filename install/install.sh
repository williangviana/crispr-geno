#!/bin/bash
# crispr-geno — One-line installer (macOS)
# Usage: curl -fsSL https://raw.githubusercontent.com/williangviana/crispr-geno/main/install/install.sh | bash

set -e

APP_NAME="crispr-geno"
REPO="williangviana/crispr-geno"
BRANCH="${CRISPR_GENO_BRANCH:-main}"
DEFAULT_DEST="$HOME/crispr-geno"

trap 'echo ""; echo "ERROR: Installation failed at the step above."; echo "Please screenshot this output and send it for help."; exit 1' ERR

echo ""
echo "============================================"
echo "  crispr-geno — Installer"
echo "============================================"
echo ""

# --- Prompt for install location (works even when piped through bash) ---
# User can press Enter for the default, or type 'b' to open a Finder picker.
DEST=""
if [ -r /dev/tty ]; then
    {
        printf "Where should crispr-geno be installed?\n"
        printf "  • Press Enter to use the default: %s\n" "$DEFAULT_DEST"
        printf "  • Type 'b' and press Enter to browse in Finder\n"
        printf "> "
    } > /dev/tty
    read USER_INPUT < /dev/tty || true

    case "$USER_INPUT" in
        b|B|browse|Browse)
            printf "Opening Finder...\n" > /dev/tty
            PARENT=$(osascript \
                -e 'try' \
                -e 'POSIX path of (choose folder with prompt "Pick the folder where crispr-geno should be installed:")' \
                -e 'on error' \
                -e 'return ""' \
                -e 'end try' 2>/dev/null || true)
            if [ -z "$PARENT" ]; then
                printf "No folder selected — using default: %s\n" "$DEFAULT_DEST" > /dev/tty
                DEST="$DEFAULT_DEST"
            else
                PARENT="${PARENT%/}"
                DEST="$PARENT/crispr-geno"
            fi
            ;;
        *)
            DEST="$DEFAULT_DEST"
            ;;
    esac
fi
DEST="${DEST:-$DEFAULT_DEST}"
DEST="${DEST/#\~/$HOME}"

echo ""
echo "Installing to: $DEST"
echo ""

# --- 1. Homebrew ---
if ! command -v brew &>/dev/null; then
    echo "[1/6] Installing Homebrew..."
    NONINTERACTIVE=1 /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
fi

if [ -x /opt/homebrew/bin/brew ]; then
    eval "$(/opt/homebrew/bin/brew shellenv)"
elif [ -x /usr/local/bin/brew ]; then
    eval "$(/usr/local/bin/brew shellenv)"
fi
echo "[1/6] Homebrew ✓"

# --- 2. Install system dependencies ---
echo "[2/6] Installing python@3.12, minimap2, samtools (may take a few minutes)..."
for pkg in python@3.12 minimap2 samtools; do
    if ! brew list "$pkg" &>/dev/null; then
        brew install "$pkg"
    fi
done
echo "[2/6] System dependencies ✓"

# Pick Python 3.12
if [ -x /opt/homebrew/opt/python@3.12/bin/python3.12 ]; then
    PY=/opt/homebrew/opt/python@3.12/bin/python3.12
elif [ -x /usr/local/opt/python@3.12/bin/python3.12 ]; then
    PY=/usr/local/opt/python@3.12/bin/python3.12
elif command -v python3.12 &>/dev/null; then
    PY=python3.12
else
    PY=python3
fi

# --- 3. Download repo to destination ---
echo "[3/6] Downloading $APP_NAME..."
if [ -e "$DEST" ] && [ -n "$(ls -A "$DEST" 2>/dev/null)" ]; then
    BACKUP="${DEST}.bak-$(date +%Y%m%d-%H%M%S)"
    echo "  Destination not empty — moving existing folder to: $BACKUP"
    mv "$DEST" "$BACKUP"
fi
mkdir -p "$DEST"
curl -fsSL "https://github.com/$REPO/archive/refs/heads/$BRANCH.tar.gz" | tar xz -C "$DEST" --strip-components=1
echo "[3/6] Downloaded ✓"

# --- 4. Virtual environment ---
cd "$DEST"
$PY -m venv .venv
# shellcheck disable=SC1091
source .venv/bin/activate
echo "[4/6] Virtual environment ✓"

# --- 5. Install Python package ---
echo "[5/6] Installing Python package..."
pip install --upgrade pip -q
pip install -e . -q
echo "[5/6] Python package ✓"

# --- 6. Link CLI on PATH ---
BIN_DIR="$HOME/.local/bin"
mkdir -p "$BIN_DIR"
ln -sf "$DEST/.venv/bin/crispr-geno" "$BIN_DIR/crispr-geno"

# Add ~/.local/bin to PATH if missing from the user's shell rc
SHELL_RC=""
case "${SHELL:-}" in
    *zsh)  SHELL_RC="$HOME/.zshrc" ;;
    *bash) SHELL_RC="$HOME/.bash_profile" ;;
    *)     SHELL_RC="$HOME/.zshrc" ;;
esac

if [ -n "$SHELL_RC" ] && ! grep -qs '\.local/bin' "$SHELL_RC"; then
    {
        echo ""
        echo "# Added by crispr-geno installer"
        echo 'export PATH="$HOME/.local/bin:$PATH"'
    } >> "$SHELL_RC"
fi
echo "[6/6] CLI linked ✓"

echo ""
echo "============================================"
echo "  Installation complete!"
echo ""
echo "  Folder:   $DEST"
echo "  Command:  ~/.local/bin/crispr-geno"
echo ""
echo "  Open a NEW Terminal window, then run:"
echo "      crispr-geno --help"
echo ""
echo "  Or double-click:"
echo "      $DEST/crispr-geno.command"
echo "============================================"
