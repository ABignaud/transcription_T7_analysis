#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import bacchus.hic as bch
import cooler
import copy
import hicstuff.hicstuff as hcs
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import ndimage
import serpentine as srp

# Import snakemake.
# mat_t7_ara30_file = snakemake.input.mat_ara30
# mat_t7_ara30_topA_file = snakemake.input.mat_ara30_topA
# mat_t7_ara60_file = snakemake.input.mat_ara60
# mat_t7_ara60_topA_file = snakemake.input.mat_ara60_topA
# cmap = snakemake.params.cmap
# res = snakemake.params.res
# res_ratio = snakemake.params.res_ratio
# width = snakemake.params.width
# cpus = snakemake.threads
mat_t7_ara30_file = '/pasteur/appa/homes/ambignau/rsg_fast/T7/data/HiC/T7_pBAD_ara30.mcool'
mat_t7_ara30_topA_file = '/pasteur/appa/homes/ambignau/rsg_fast/T7/data/HiC/T7_pBAD_ara30_topA.mcool'
mat_t7_ara60_file = '/pasteur/appa/homes/ambignau/rsg_fast/T7/data/HiC/T7_pBAD_ara60.mcool'
mat_t7_ara60_topA_file = '/pasteur/appa/homes/ambignau/rsg_fast/T7/data/HiC/T7_pBAD_ara60_topA.mcool'
cmap = 'Reds'
res = 1000
res_ratio = 5000
width = 64000
cpus = 25

# outfile = str(snakemake.output.matrices)
# out_zoom = str(snakemake.output.zoom)
# out_zoom_ratio = str(snakemake.output.zoom_ratio)

outfile = str('/pasteur/appa/homes/ambignau/rsg_fast/T7/data/figures/T7_system/TopA/1.pdf')
out_zoom = str('/pasteur/appa/homes/ambignau/rsg_fast/T7/data/figures/T7_system/TopA/2.pdf')
out_zoom_ratio = str('/pasteur/appa/homes/ambignau/rsg_fast/T7/data/figures/T7_system/TopA/3.pdf')

# Make sure output directory exists.
os.makedirs(str('/pasteur/appa/homes/ambignau/rsg_fast/T7/data/figures/T7_system/TopA/'), exist_ok=True)

# Import matrices
mat_t7_ara30 = cooler.Cooler(f"{mat_t7_ara30_file}::/resolutions/{res}").matrix(
    balance=True, sparse=False
)[:]
mat_t7_ara30[np.isnan(mat_t7_ara30)] = 0
mat_t7_ara30_topA = cooler.Cooler(f"{mat_t7_ara30_topA_file}::/resolutions/{res}").matrix(
    balance=True, sparse=False
)[:]
mat_t7_ara30_topA[np.isnan(mat_t7_ara30_topA)] = 0
mat_t7_ara60 = cooler.Cooler(f"{mat_t7_ara60_file}::/resolutions/{res}").matrix(
    balance=True, sparse=False
)[:]
mat_t7_ara60[np.isnan(mat_t7_ara60)] = 0
mat_t7_ara60_topA = cooler.Cooler(f"{mat_t7_ara60_topA_file}::/resolutions/{res}").matrix(
    balance=True, sparse=False
)[:]
mat_t7_ara60_topA[np.isnan(mat_t7_ara60_topA)] = 0

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

mat_t7_ara30_sparse = cooler.Cooler(
    f"{mat_t7_ara30_file}::/resolutions/{res_ratio}"
).matrix(
    balance=False, sparse=True
)[:]
mat_t7_ara30_topA_sparse = cooler.Cooler(
    f"{mat_t7_ara30_topA_file}::/resolutions/{res_ratio}"
).matrix(
    balance=False, sparse=True
)[:]
mat_t7_ara60_sparse = cooler.Cooler(
    f"{mat_t7_ara60_file}::/resolutions/{res_ratio}"
).matrix(
    balance=False, sparse=True
)[:]
mat_t7_ara60_topA_sparse = cooler.Cooler(
    f"{mat_t7_ara60_topA_file}::/resolutions/{res_ratio}"
).matrix(
    balance=False, sparse=True
)[:]

(
    _,
    _,
    ratio_log_ara30,
    ratio_srp_ara30,
    subsample,
) = contact_map_ratio(
    mat_t7_ara30_sparse,
    mat_t7_ara30_topA_sparse,
    iterations=25,
    threshold=100,
    cpus=cpus,
)
(
    _,
    _,
    ratio_log_ara60,
    ratio_srp_ara60,
    subsample,
) = contact_map_ratio(
    mat_t7_ara60_sparse,
    mat_t7_ara60_topA_sparse,
    iterations=25,
    threshold=100,
    cpus=cpus,
)

###
# PLOT MATRICES
###

title = ["Ara 30 min", "Ara 60 min ", "Ara 30 min - topA", "Ara 60 min - topA"]
hic = [mat_t7_ara30, mat_t7_ara60, mat_t7_ara30_topA, mat_t7_ara60_topA]
ratio = [
    ratio_srp_ara30 - ratio_srp_ara30.mean(),
    ratio_srp_ara60 - ratio_srp_ara60.mean(),
]
start = 0
end = len(hic[0])

fig, ax = plt.subplots(2, 3, figsize=(10, 10))
for i in range(2):
    for j in range(2):
        ax[i, j].imshow(
            hic[i + 2 * j][start:end, start:end],
            vmin=0,
            vmax=np.nanpercentile(hic[0], 98),
            cmap="Reds",
            extent=(start, end, end, start),
        )
        ax[i, j].set_title(title[i + 2 * j])
    ax[i, 2].imshow(
        ratio[i][start:end//5, start:end//5],
        vmin=-2,
        vmax=2,
        cmap="seismic",
        extent=(start, end, end, start),
    )
    ax[i, 2].set_title("topA+/topA-")

plt.savefig(outfile)
plt.close()


###
# Plot zoom without RNAseq.
###

# Rotate matrices.
hic_rot = copy.copy(hic)
for i in range(len(hic)):
    hic_rot[i] = bch.interpolate_white_lines(hic_rot[i])
    hic_rot[i][np.isnan(hic_rot[i])] = 0
    hic_rot[i] = ndimage.rotate(hic_rot[i], 45, reshape=True)

fig, ax = plt.subplots(
    1, 4, figsize=(20, 10)
)

start = 300_000
end = 480_000

row1 = len(hic_rot[0]) // 2 - int(np.sqrt(2) * (width / res))
row2 = len(hic_rot[0]) // 2 + int(np.sqrt(2) * (width / res))
col1 = int((start // res) * np.sqrt(2))
col2 = int((end // res) * np.sqrt(2))

for i in range(4):
    # Hicmap
    im = ax[i].imshow(
        hic_rot[i][row1:row2, col1:col2],
        cmap="Reds",
        vmin=0,
        vmax=np.nanpercentile(hic[0], 99),
        extent=(col1 / np.sqrt(2), col2 / np.sqrt(2), width // res, -width // res),
    )
    ax[i].get_xaxis().set_visible(False)
    ax[i].set_title(title[i], size=18)
    if i == 0:
        ax[i].tick_params(axis="both", labelsize=16)
        ax[i].set_ylabel("Genomic coordinates (kb)", size=16)
    else:
        ax[i].get_yaxis().set_visible(False)

# Colorbar
cbar = fig.colorbar(im, ax=ax.ravel().tolist(), shrink=0.5, anchor=(1.2, 0.8))
cbar.ax.tick_params(labelsize=16)

# Adjust space between plot
plt.subplots_adjust(wspace=0.4, hspace=0.0, left=0.05, right=0.88)
plt.savefig(out_zoom)
plt.close()

###
# Plot zoom ratio.
###

start, end = 0 // 5, 750 // 5
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
