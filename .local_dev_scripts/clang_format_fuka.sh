# This script formats the source code related to FUKA using clang-format
# 1. No changes to Kadath source code is made
# 2. EOS Standalone codes that are from different repositories
#    are also not formatted
# 3. Explicit file paths are used to avoid formatting unrelated files

#!/usr/bin/env bash
set -euo pipefail

source $HOME_KADATH/.venv/bin/activate

# Directories to format recursively
DIRS=(
  codes/FUKA
  codes/PythonTools
  include/Configurator
  include/EOS/ghl_eos_helpers/
  include/FUKA_Solvers
  include/Solvers
  src/Utilities/Configurator
  src/Utilities/EOS
  src/Utilities/Exporters
  src/Utilities/PN
)

# File extensions to format
find "${DIRS[@]}" -type f \( \
  -name "*.c" -o -name "*.cc" -o -name "*.cpp" -o \
  -name "*.h" -o -name "*.hh" -o -name "*.hpp" \
\) -print0 | xargs -0 -r clang-format -i

# Re-stage formatted changes
git add "${DIRS[@]}"


FILES=(
  include/bco_utilities.hpp
  include/coord_fields.hpp
  include/exporter_utilities.hpp
  include/EOS/EOS.hh
  include/EOS/FUKA_EOS_Utilities.hh
  include/EOS/FUKA_EOS_Wrapper.hh
  include/EOS/standalone/EOS_parfile_parser.hpp
  include/EOS/standalone/tov.hh
)
for file in "${FILES[@]}"; do
  clang-format -i "$file"
done

git add "${FILES[@]}"