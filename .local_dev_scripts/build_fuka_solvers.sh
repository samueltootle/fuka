#!/usr/bin/env bash
set -euo pipefail

# File extensions to format
COMPILE_PATHS=$(find $HOME_KADATH/codes/FUKA -name "compile")

for f in ${COMPILE_PATHS[@]}; do
  echo "Compiling $f"
  cd "$(dirname "$f")"
  bash "$(basename "$f")"
  cd -
done