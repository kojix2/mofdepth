#!/usr/bin/env bash
set -euo pipefail

# This script regenerates mosdepth fixtures used by specs.
# Requirements:
#   - MOSDEPTH_PATH env var pointing to mosdepth binary (or edit path below)
#   - Test inputs located under ./mosdepth/tests as in repo

MOSDEPTH_BIN=${MOSDEPTH_PATH:-"$(pwd)/mosdepth/mosdepth"}
FIXDIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$FIXDIR/.." && pwd)"
TMPDIR="${TMPDIR:-/tmp/moffdepth_fixtures}"

mkdir -p "$TMPDIR"

# Inputs
BAM_OVL="$ROOT/../mosdepth/tests/ovl.bam"
BAM_FRAG="$ROOT/../mosdepth/tests/full-fragment-pairs.bam"
BED_TRACK="$ROOT/../mosdepth/tests/track.bed"

if [[ ! -x "$MOSDEPTH_BIN" ]]; then
  echo "mosdepth binary not found: $MOSDEPTH_BIN" >&2
  exit 1
fi

run_case() {
  local prefix="$1"; shift
  local bam="$1"; shift
  local args=("$@")
  echo "Generating: $prefix with args: ${args[*]}"
  (cd "$TMPDIR" && "$MOSDEPTH_BIN" "$prefix" "${args[@]}" "$bam")
  # Move outputs to fixtures, normalize names to prefix basename
  mv "$TMPDIR/$(basename "$prefix")."* "$FIXDIR/"
}

# Basic per-base + summary/dist
run_case "ovl" "$BAM_OVL"

# Flags
run_case "ovl_flag" "$BAM_OVL" -F 4
run_case "ovl_include" "$BAM_OVL" -i 2

# Quantize
run_case "ovl_quant" "$BAM_OVL" -q 0:1:1000

# Regions (window and bed)
run_case "ovl_win" "$BAM_OVL" -b 100
if [[ -f "$BED_TRACK" ]]; then
  run_case "ovl_bed" "$BAM_OVL" -b "$BED_TRACK"
fi

# Thresholds
run_case "ovl_thresh" "$BAM_OVL" -T 0,1,2 -b 100

# Fast
run_case "ovl_fast" "$BAM_OVL" -x

# Fragment
if [[ -f "$BAM_FRAG" ]]; then
  run_case "ovl_frag" "$BAM_FRAG" -a
fi

echo "Fixtures written to $FIXDIR"
