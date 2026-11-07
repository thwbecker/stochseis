#!/bin/bash
# Copyright (C) 2026 Thorsten Becker. Written by Claude (Anthropic) with Thorsten Becker.
# SPDX-License-Identifier: GPL-3.0-or-later
# Reproduce the single-fault experiments of docs/experiments.md.
# usage: run_all.sh [nproc] [outdir]
nproc=${1-4}
outdir=${2-experiments_out}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$outdir"
for spec in "$here"/e*.json
do
    name=$(basename "$spec" .json)
    python -m stochseis.cli sweep "$spec" "$outdir/$name" --nproc "$nproc" --no-catalogs
done
