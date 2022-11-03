#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pickletools import markobject
import bacchus.hic as bch
import bacchus.io as bcio
import cooler
import copy
import hicstuff.hicstuff as hcs
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
import serpentine as srp

# Import snakemake.
mat_t7_glu = snakemake.input.mat_t7_glu
mat_t7_ara = snakemake.input.mat_t7_ara
mat_t7_glu_rif = snakemake.input.mat_t7_glu_rif
mat_t7_ara_rif = snakemake.input.mat_t7_ara_rif
rna_t7_glu = snakemake.input.rna_t7_glu
rna_t7_ara = snakemake.input.rna_t7_ara
rna_t7_glu_rif = snakemake.input.rna_t7_glu_rif
rna_t7_ara_rif = snakemake.input.rna_t7_ara_rif

cmap = snakemake.params.cmap
res = snakemake.params.res
width = snakemake.params.width

outfile = str(snakemake.output.matrices)
out_zoom = str(snakemake.output.zoom)
out_zoom_ratio = str(snakemake.output.zoom_ratio)

###
# PLOT MATRICES
###

# Import matrices
mat_t7_glu = cooler.Cooler(f'{mat_t7_glu}::/resolutions/{res}').matrix(balance=False, sparse=True)[:]
mat_t7_ara = cooler.Cooler(f'{mat_t7_ara}::/resolutions/{res}').matrix(balance=False, sparse=True)[:]
mat_t7_glu_rif = cooler.Cooler(f'{mat_t7_glu_rif}::/resolutions/{res}').matrix(balance=False, sparse=True)[:]
mat_t7_ara_rif = cooler.Cooler(f'{mat_t7_ara_rif}::/resolutions/{res}').matrix(balance=False, sparse=True)[:]

# Do the ratio
def contact_map_ratio(
    M1,
    M2,
    iterations=10,
    threshold=10,
    cpus=16,
):
    """Function to compute the log ratio of two matrices and use serpentine to
    return a smoother log ratio.
    """
    M1 = M1.tocoo()
    M2 = M2.tocoo()
    m1 = M1.sum()
    m2 = M2.sum()
    if m1 < m2:
        M2 = hcs.subsample_contacts(M2, int(m1))
    else:
        M1 = hcs.subsample_contacts(M1, int(m2))
    subsample = min(m1, m2)
    print("Matrices have been subsampled to {0} contacts.".format(subsample))
    # M1 = M1.tocsr() - M1.tocsr().multiply(
    #     get_win_density(M1, win_size=3, sym_upper=True) < 2 / 9
    # )
    # M2 = M2.tocsr() - M2.tocsr().multiply(
    #     get_win_density(M2, win_size=3, sym_upper=True) < 2 / 9
    # )
    M1 = bch.get_symmetric(M1).toarray()
    M2 = bch.get_symmetric(M2).toarray()
    if threshold == "auto":
        trend, threshold = srp.MDbefore(M1, M2, show=False)
    else:
        threshold = int(threshold)
    M1_serp, M2_serp, ratio_srp = srp.serpentin_binning(
        M1,
        M2,
        parallel=cpus,
        triangular=False,
        threshold=threshold,
        minthreshold=threshold / 5,
        iterations=iterations,
        verbose=False,
    )
    ratio_log = np.log2(M2) - np.log2(M1)
    return M1_serp, M2_serp, ratio_log, ratio_srp, subsample

mat_t7_glu_rif, mat_t7_glu, ratio_log_glu, ratio_srp_glu, subsample = contact_map_ratio(mat_t7_glu_rif, mat_t7_glu)
mat_t7_ara_rif, mat_t7_ara, ratio_log_ara, ratio_srp_ara, subsample = contact_map_ratio(mat_t7_ara_rif, mat_t7_ara)

# Plot it
title = ["Glu", "Glu Rif", "Ara", "Ara Rif"]
hic = [mat_t7_glu, mat_t7_glu_rif, mat_t7_ara, mat_t7_ara_rif]
ratio = [ratio_srp_glu - ratio_srp_glu.mean(), ratio_srp_ara - ratio_srp_ara.mean()]
start = 0
end = len(hic[0])

fig, ax = plt.subplots(2,3, figsize=(10,10))
for i in range(2):
    for j in range(2):
        ax[i,j].imshow(
            hic[i+2*j][start:end, start:end],
            vmin=0,
            vmax= np.nanpercentile(hic[0], 99),
            cmap="Reds",
            extent=(start, end ,end, start),
        )
        ax[i,j].set_title(title[i+2*j])
    ax[i,2].imshow(
        ratio[i][start:end, start:end],
        vmin=-2,
        vmax=2,
        cmap="seismic",
        extent=(start, end ,end, start),
    )
    ax[i,2].set_title('Rif+/Rif-')

plt.savefig(outfile)

###
# Plot zoom.
###

# Rotate matrices.
hic_rot = copy.copy(hic)
for i, mat in enumerate(hic):
    hic_rot[i] = bch.interpolate_white_lines(hic[i])
    hic_rot[i] = ndimage.rotate(hic_rot[i], 45, reshape=True)

# Import RNA
rna_binning = 100
rna_t7_glu, _ = bcio.extract_big_wig(rna_t7_glu, rna_binning)
rna_t7_ara, _ = bcio.extract_big_wig(rna_t7_ara, rna_binning)
rna_t7_glu_rif, _ = bcio.extract_big_wig(rna_t7_glu_rif, rna_binning)
rna_t7_ara_rif, _ = bcio.extract_big_wig(rna_t7_ara_rif, rna_binning)
rna = [rna_t7_glu, rna_t7_glu_rif, rna_t7_ara, rna_t7_ara_rif]

fig, ax = plt.subplots(2, 4, figsize=(10,10), gridspec_kw={'height_ratios': [8, 4]}, sharex=False)

start = 300_000
end = 480_000

row1 = len(hic_rot[0]) // 2 - int(np.sqrt(2) * width)
row2 = len(hic_rot[0]) // 2 + int(np.sqrt(2) * width)
col1 = int((start // res) * np.sqrt(2))
col2 = int((end // res) * np.sqrt(2))

for i in range(4):
    # Hicmap
    im = ax[0, i].imshow(
        hic_rot[i][row1:row2, col1:col2],
        cmap="Reds",
        vmax=np.percentile(hic[0], 99.9),
        extent=(col1 / np.sqrt(2), col2 / np.sqrt(2), width, -width),
    )
    ax[0, i].get_xaxis().set_visible(False)
    ax[0, i].tick_params(axis='both', labelsize=16)
    ax[0, i].set_ylabel("Genomic coordinates (kb)", size=16)
    
    # RNAseq
    ax[1, i].set_xlabel("Genomic coordinates (kb)", size=16)
    ax[1, i].set_ylabel("Transcription (CPM)", size=16)
    ax_rna = ax[1, i]
    p1 = ax_rna.plot(
        np.arange(start, end, 1 / rna_binning),
        rna[i][int(start * rna_binning):int(end * rna_binning)],
        color='k',
        label="RNAseq",
        alpha=.7,
        linewidth=1,
    )
    ax_rna.set_ylim(0, 8000)
    ax_rna.tick_params(axis='both', labelsize=16, color='k', labelcolor='k')
    ax_rna.spines['top'].set_visible(False)
    ax_rna.spines['right'].set_visible(False)
    
# Colorbar
cbar = fig.colorbar(im, ax=ax.ravel().tolist(), shrink=.5, anchor=(1.2, .8))
cbar.ax.tick_params(labelsize=16)

# Adjust space between plot
plt.subplots_adjust(wspace=.4, hspace=.0, left=0.05, right=0.88)
plt.savefig(out_zoom)

###
# Plot zoom ratio.
###

