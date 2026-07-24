#!/usr/bin/env bash
set -euo pipefail

lock_dir="${R_LIBS_USER:-}/00LOCK-coord.serial"
if [[ -n "${R_LIBS_USER:-}" && -d "$lock_dir" ]]; then
	rm -rf "$lock_dir"
fi

R -q -e "devtools::install('.')"

echo "coord.serial installed. Start R and run: library(coord.serial)"
