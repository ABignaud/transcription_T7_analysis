#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Script to plot the matrix from a cool file using bacchus functions.

import cooler
import bacchus.plot as bcp
import os

# Create outdir if necessary.
os.makedirs(str(snakemake.params.outdir), exist_ok=True)

# Import the matrix.
mat = cooler.Cooler(
    f"{str(snakemake.input.mat)}::/resolutions/{snakemake.params.res}"
).matrix(balance=True, sparse=False)[:]

# Plot the matrix.
bcp.contact_map(
    mat,
    axis="kb",
    binning=snakemake.params.res,
    cmap=snakemake.params.cmap,
    out_file=str(snakemake.output),
    title=snakemake.params.title,
)
