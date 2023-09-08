#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import bacchus.io as bcio
import bacchus.plot as bcp
import bacchus.hic as bch
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import os
from os.path import dirname

# Make sure that the figures directory exists.
os.makedirs(dirname(snakemake.output.mat_fig), exist_ok=True)
os.makedirs(dirname(snakemake.output.corr_fig), exist_ok=True)

# Import parameters
binning = int(snakemake.params.binning)
label = str(snakemake.params.species)
circular = str(snakemake.params.circular)

# Import contact map as a dense matrix.
M = bcio.build_map(
    matrix_files=[str(snakemake.input.hic)],
    fragment_file="none",
    bin_size=binning,
    mat_format="cool",
    normalize=True,
    subsample=0,
)

# Plot the contact map.
bcp.contact_map(
    M,
    axis="Mb",
    binning=binning,
    cmap="Reds",
    dpi=500,
    out_file=str(snakemake.output.mat_fig),
    title=label,
    vmax=99,
)

# Import RNA tracks.
rna, chrom_start = bcio.extract_big_wig(
    file=str(snakemake.input.rna),
    binning=binning,
    circular=circular,
    sigma=None,
    ztransform=None,
)
# Ugly loop to have the chrom length...
chrom_start_size = {}
for i, name in enumerate(chrom_start):
    if i != 0:
        chrom_start_size[prev_name] = {
            "start": start,
            "length": chrom_start[name] - start,
        }
    prev_name = name
    start = chrom_start[name]
chrom_start_size[prev_name] = {
    "start": start,
    "length": chrom_start[name] - start,
}

# Compute HiC signal.
hic = bch.compute_hic_signal(M, binning=binning, start=5000, stop=10000)
print(st.spearmanr(hic, rna[:-2])[0])

# Plot correlation between HiC signal and RNAseq.
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
ax.scatter(hic, np.log(rna[:-2]))
ax.set_xlabel("HiC signal")
ax.set_ylabel("Transcription (log)")
ax.set_title(label)
ax.text(
    x=np.nanpercentile(hic, 99),
    y=np.nanpercentile(np.log(rna), 1),
    s=f"corr={st.spearmanr(hic, rna[:-2])[0]:.2f}",
)
plt.savefig(snakemake.output.corr_fig, dpi=100)
