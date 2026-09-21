#!/bin/sh
# Runs every Stage 3 deck in its own folder.
#   ./run_all.sh /path/to/ls-dyna 4
set -e
SOLVER="${1:-ls-dyna}"
NCPU="${2:-4}"
HERE="$(cd "$(dirname "$0")" && pwd)"

echo "=== stage3_lance_only ==="
(cd "$HERE/stage3_lance_only" && "$SOLVER" i=stage3_lance_only.k ncpu="$NCPU" memory=200m)

echo "=== stage3_extract_mu000 ==="
(cd "$HERE/stage3_extract_mu000" && "$SOLVER" i=stage3_extract_mu000.k ncpu="$NCPU" memory=200m)

echo "=== stage3_extract_mu020 ==="
(cd "$HERE/stage3_extract_mu020" && "$SOLVER" i=stage3_extract_mu020.k ncpu="$NCPU" memory=200m)

echo "=== stage3_extract_mu030 ==="
(cd "$HERE/stage3_extract_mu030" && "$SOLVER" i=stage3_extract_mu030.k ncpu="$NCPU" memory=200m)

echo "=== stage3_insert_mu020 ==="
(cd "$HERE/stage3_insert_mu020" && "$SOLVER" i=stage3_insert_mu020.k ncpu="$NCPU" memory=200m)

echo "=== stage3_extract_tpa_g010 ==="
(cd "$HERE/stage3_extract_tpa_g010" && "$SOLVER" i=stage3_extract_tpa_g010.k ncpu="$NCPU" memory=200m)

echo "=== stage3_extract_tpa_g085 ==="
(cd "$HERE/stage3_extract_tpa_g085" && "$SOLVER" i=stage3_extract_tpa_g085.k ncpu="$NCPU" memory=200m)

echo "=== stage3_extract_mu020_fine ==="
(cd "$HERE/stage3_extract_mu020_fine" && "$SOLVER" i=stage3_extract_mu020_fine.k ncpu="$NCPU" memory=200m)

echo "=== stage3_extract_mu020_slow ==="
(cd "$HERE/stage3_extract_mu020_slow" && "$SOLVER" i=stage3_extract_mu020_slow.k ncpu="$NCPU" memory=200m)

echo 'Done. Next:'
echo '  python scripts/extract_stage3.py "'"$HERE"'"'
echo '  python scripts/postprocess_stage3.py'
