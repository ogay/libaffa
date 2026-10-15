#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
binary=$(mktemp "${TMPDIR:-/tmp}/libaffa-abs.XXXXXX")
trap 'rm -f "$binary"' EXIT HUP INT TERM
"${CXX:-c++}" -std=c++11 -O2 -Isrc src/*.cpp tests/regression_abs.cpp -o "$binary"
"$binary"
