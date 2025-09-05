This directory holds pre-generated reference outputs produced by mosdepth for test BAMs used in specs.

How to refresh fixtures:
- Use `mosdepth.sh` in this directory. It will run mosdepth against known inputs and store outputs here.
- Ensure `MOSDEPTH_PATH` points to a working mosdepth binary (or adjust the script).

Files naming convention:
- <name>.<suffix> as produced by mosdepth, e.g. `ovl.per-base.bed.gz`, `ovl.mosdepth.summary.txt`, etc.
