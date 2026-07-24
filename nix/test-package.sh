#!/usr/bin/env bash
set -euo pipefail

# Build, install, and run tests for the local package source.
PKG_NAME="coord.serial"
PKG_VERSION="$(awk -F': *' '$1=="Version" {print $2}' DESCRIPTION)"
TARBALL="${PKG_NAME}_${PKG_VERSION}.tar.gz"
LOCAL_LIB="$(mktemp -d "${TMPDIR:-/tmp}/coord-serial-lib.XXXXXX")"
trap 'rm -rf "$LOCAL_LIB"' EXIT

export R_LIBS_USER="$LOCAL_LIB"

rm -f "$TARBALL"
R CMD build .
R CMD INSTALL -l "$LOCAL_LIB" "$TARBALL"
R -q -e ".libPaths(c('${LOCAL_LIB}', .libPaths())); testthat::test_dir('tests/testthat', load_package = 'source')"

echo "Completed: built, installed, and tested ${TARBALL}"
