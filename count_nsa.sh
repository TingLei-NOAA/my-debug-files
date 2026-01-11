  #!/bin/bash
  # usage: ./count_nsa.sh "DataFix/bump_401km_0p04sigma_mpi160/bumploc_401km0p04sigma*.nc"

  set -euo pipefail

  if [ $# -ne 1 ]; then
    echo "Usage: $0 \"pattern\""
    exit 1
  fi

  pattern=$1

  if ! ls -1 $pattern >/dev/null 2>&1; then
    echo "No files match pattern: $pattern"
    exit 1
  fi

  grand_total=0

  while IFS= read -r f; do
    mapfile -t nsa_lines < <(ncdump -h "$f" | grep -F "nsa =")
    count=${#nsa_lines[@]}

    if [ "$count" -ne 2 ]; then
      echo "ERROR: $f has $count occurrences of \"nsa =\" (expected 2)"
      exit 2
    fi

    second_line=${nsa_lines[1]}
    value=$(awk -F= '{val=$2; gsub(/[ ;]/,"",val); print val}' <<<"$second_line")

    echo "$f: $value"
    grand_total=$((grand_total + value))
  done < <(ls -1 $pattern 2>/dev/null)

  echo "Grand total: $grand_total"
