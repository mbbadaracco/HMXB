#!/usr/bin/env bash
# Install a user-space TeX Live into work/texlive.
#
# Why: this machine has no system LaTeX (the conda environment
# ~/Software/gwtex ships TeX binaries but zero packages -- its
# share/texmf-dist/tex/latex is empty), and sudo requires a password, so
# apt is not an option.  This script is the sequence that was actually
# run; re-running it reproduces the toolchain the summary is built
# with.
set -euo pipefail

W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

curl -sSL -o "$TMP/install-tl-unx.tar.gz" \
  https://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz
tar xzf "$TMP/install-tl-unx.tar.gz" -C "$TMP"
INST="$(echo "$TMP"/install-tl-*/install-tl)"

cat > "$TMP/tl.profile" <<PROFILE
selected_scheme scheme-small
TEXDIR $W/texlive
TEXMFLOCAL $W/texlive/texmf-local
TEXMFSYSVAR $W/texlive/texmf-var
TEXMFSYSCONFIG $W/texlive/texmf-config
TEXMFVAR $W/texlive/user-texmf-var
TEXMFCONFIG $W/texlive/user-texmf-config
TEXMFHOME $W/texlive/texmf-home
instopt_adjustpath 0
instopt_adjustrepo 1
instopt_letter 0
instopt_portable 0
instopt_write18_restricted 1
tlpdbopt_autobackup 0
tlpdbopt_install_docfiles 0
tlpdbopt_install_srcfiles 0
PROFILE

"$INST" -profile "$TMP/tl.profile"

# Packages added on top of scheme-small.
export PATH="$W/texlive/bin/x86_64-linux:$PATH"
tlmgr install sttools txfonts collection-latexextra \
              collection-fontsrecommended collection-bibtexextra

echo "TeX Live installed at $W/texlive"
"$W/texlive/bin/x86_64-linux/pdflatex" --version | head -1
