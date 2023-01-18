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
import os
from scipy import ndimage
import serpentine as srp

# Import snakemake.
mat_t7_glu_file = snakemake.input.mat_t7_glu
mat_t7_ara_file = snakemake.input.mat_t7_ara
mat_t7_glu_rif_file = snakemake.input.mat_t7_glu_rif
mat_t7_ara_rif_file = snakemake.input.mat_t7_ara_rif
rna_t7_glu = snakemake.input.rna_t7_glu
rna_t7_ara = snakemake.input.rna_t7_ara
rna_t7_glu_rif = snakemake.input.rna_t7_glu_rif
rna_t7_ara_rif = snakemake.input.rna_t7_ara_rif
chip_t7_glu = snakemake.input.chip_t7_glu
chip_t7_ara = snakemake.input.chip_t7_ara
chip_t7_glu_rif = snakemake.input.chip_t7_glu_rif
chip_t7_ara_rif_R1 = snakemake.input.chip_t7_ara_rif_R1
chip_t7_ara_rif_R2 = snakemake.input.chip_t7_ara_rif_R2

cmap = snakemake.params.cmap
res = snakemake.params.res
res_ratio = snakemake.params.res_ratio
width = snakemake.params.width
cpus = snakemake.threads

outdir = str(snakemake.params.outdir)
outfile = str(snakemake.output.matrices)
out_zoom = str(snakemake.output.zoom)
out_zoom_ratio = str(snakemake.output.zoom_ratio)

# Make sure output directory exists.
os.makedirs(outdir, exist_ok=True)

# Import matrices
mat_t7_glu = cooler.Cooler(f"{mat_t7_glu_file}::/resolutions/{res}").matrix(
    balance=True, sparse=False
)[:]
mat_t7_glu[np.isnan(mat_t7_glu)] = 0
mat_t7_ara = cooler.Cooler(f"{mat_t7_ara_file}::/resolutions/{res}").matrix(
    balance=True, sparse=False
)[:]
mat_t7_ara[np.isnan(mat_t7_ara)] = 0
mat_t7_glu_rif = cooler.Cooler(f"{mat_t7_glu_rif_file}::/resolutions/{res}").matrix(
    balance=True, sparse=False
)[:]
mat_t7_glu_rif[np.isnan(mat_t7_glu_rif)] = 0
mat_t7_ara_rif = cooler.Cooler(f"{mat_t7_ara_rif_file}::/resolutions/{res}").matrix(
    balance=True, sparse=False
)[:]
mat_t7_ara_rif[np.isnan(mat_t7_ara_rif)] = 0

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

mat_t7_glu_sparse = cooler.Cooler(f"{mat_t7_glu_file}::/resolutions/{res_ratio}").matrix(
    balance=False, sparse=True
)[:]
mat_t7_ara_sparse = cooler.Cooler(f"{mat_t7_ara_file}::/resolutions/{res_ratio}").matrix(
    balance=False, sparse=True
)[:]
mat_t7_glu_rif_sparse = cooler.Cooler(f"{mat_t7_glu_rif_file}::/resolutions/{res_ratio}").matrix(
    balance=False, sparse=True
)[:]
mat_t7_ara_rif_sparse = cooler.Cooler(f"{mat_t7_ara_rif_file}::/resolutions/{res_ratio}").matrix(
    balance=False, sparse=True
)[:]

(
    _,
    _,
    ratio_log_glu,
    ratio_srp_glu,
    subsample,
) = contact_map_ratio(mat_t7_glu_sparse, mat_t7_ara_sparse, iterations=25, threshold=100, cpus=cpus)
(
    _,
    _,
    ratio_log_ara,
    ratio_srp_ara,
    subsample,
) = contact_map_ratio(mat_t7_glu_rif_sparse, mat_t7_ara_rif_sparse, iterations=25, threshold=100, cpus=cpus)


###
# PLOT MATRICES
###

title = ["Glucose", "Arabinose", "Glucose + Rif", "Arabinose + Rif"]
hic = [mat_t7_glu, mat_t7_ara, mat_t7_glu_rif, mat_t7_ara_rif]
ratio = [
    ratio_srp_glu - ratio_srp_glu.mean(),
    ratio_srp_ara - ratio_srp_ara.mean(),
]
start = 0
end = len(hic[0])

fig, ax = plt.subplots(2, 3, figsize=(10, 10))
for i in range(2):
    for j in range(2):
        ax[i, j].imshow(
            hic[i + 2 * j][start:end, start:end],
            vmin=0,
            vmax=np.nanpercentile(hic[0], 99),
            cmap="Reds",
            extent=(start, end, end, start),
        )
        ax[i, j].set_title(title[i + 2 * j])
    ax[i, 2].imshow(
        ratio[i][start:end, start:end],
        vmin=-2,
        vmax=2,
        cmap="seismic",
        extent=(start, end, end, start),
    )
    ax[i, 2].set_title("Rif+/Rif-")

plt.savefig(outfile, dpi=300)
plt.close()

###
# Plot zoom ratio.
###

start, end = 260 // 5, 560 // 5
fig, ax = plt.subplots(2, 1, figsize=(5, 10))
for i in range(2):
    ax[i].imshow(
        ratio[i][start:end, start:end],
        vmin=-2,
        vmax=2,
        cmap="seismic",
        extent=(start, end, end, start),
    )
plt.savefig(out_zoom_ratio)
plt.close()

###
# Plot zoom.
###

# Rotate matrices.
hic_rot = copy.copy(hic)
for i in range(len(hic)):
    hic_rot[i] = bch.interpolate_white_lines(hic_rot[i])
    hic_rot[i][np.isnan(hic_rot[i])] = 0
    hic_rot[i] = ndimage.rotate(hic_rot[i], 45, reshape=True)

start = 240_000
end = 420_000

row1 = len(hic_rot[0]) // 2 - int(np.sqrt(2) * (width / res))
row2 = len(hic_rot[0]) // 2 + int(np.sqrt(2) * (width / res))
col1 = int((start // res) * np.sqrt(2))
col2 = int((end // res) * np.sqrt(2))

# Import RNA
rna_binning = 100
rna_t7_glu, _ = bcio.extract_big_wig(rna_t7_glu, rna_binning)
rna_t7_ara, _ = bcio.extract_big_wig(rna_t7_ara, rna_binning)
rna_t7_glu_rif, _ = bcio.extract_big_wig(rna_t7_glu_rif, rna_binning)
rna_t7_ara_rif, _ = bcio.extract_big_wig(rna_t7_ara_rif, rna_binning)
rna = [rna_t7_glu, rna_t7_glu_rif, rna_t7_ara, rna_t7_ara_rif]

# Import ChIP
chip_binning = 100
chip_t7_glu, _ = bcio.extract_big_wig(chip_t7_glu, chip_binning)
chip_t7_ara, _ = bcio.extract_big_wig(chip_t7_ara, chip_binning)
chip_t7_glu_rif, _ = bcio.extract_big_wig(chip_t7_glu_rif, chip_binning)
chip_t7_ara_rif_R1, _ = bcio.extract_big_wig(chip_t7_ara_rif_R1, chip_binning)
chip_t7_ara_rif_R2, _ = bcio.extract_big_wig(chip_t7_ara_rif_R2, chip_binning)
T7_pol = [chip_t7_glu, chip_t7_glu_rif, chip_t7_ara, chip_t7_ara_rif_R1]

fig, ax = plt.subplots(
    3,
    4,
    figsize=(20, 6),
    gridspec_kw={"height_ratios": [15, 7, 7]},
    sharex=True,
)
for j in range(4):

    # HiC
    im = ax[0, j].imshow(
        hic_rot[j][row1:row2, col1:col2],
        cmap="Reds",
        vmax=0.002,
        vmin=0,
        extent=(col1 / np.sqrt(2), col2 / np.sqrt(2), width // res, -width // res),
    )

    # RNAseq
    ax[1, j].fill_between(
        np.arange(start / res, end / res, 0.1), rna[j][start // 100 : end // 100], color="k"
    )
    ax[1, j].set_ylim(0, 5000)

    # T7 pol Chip_seq
    ax[2, j].plot(
        np.arange(start / res, end / res, 0.1),
        T7_pol[j][start  // 100:end // 100],
    )
    ax[2, j].set_ylim(0, 5000)

    # Legend
    ax[0, j].set_title(title[j], size=16)
    ax[0, j].tick_params(axis="both", labelsize=12)
    ax[2, j].tick_params(axis="both", labelsize=12)
    if j != 0:
        ax[0, j].get_yaxis().set_visible(False)

    ax[1, j].set_xlabel("Genomic coordinates (kb)", fontsize=14)
ax[0, 0].set_ylabel("Genomic\ndistance (kb)", fontsize=14)
ax[1, 0].set_ylabel("Transcription\n(CPM)", fontsize=14)


# Colorbar
cbar = fig.colorbar(im, ax=ax.ravel().tolist(), shrink=0.3, anchor=(-0.2, 0.4))
cbar.ax.tick_params(labelsize=15)
plt.savefig(out_zoom)
