#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import bacchus.hic as bch
import cooler
import hicstuff.hicstuff as hcs
import matplotlib.pyplot as plt
import numpy as np
import copy 
import os
from scipy import ndimage
import serpentine as srp

# Import snakemake.
mat_t7_ara_file = snakemake.input.mat_ara
mat_t7_novo_file = snakemake.input.mat_novo
mat_t7_ara_rif_file = snakemake.input.mat_ara_rif
mat_t7_novo_rif_file = snakemake.input.mat_novo_rif
cmap = snakemake.params.cmap
res = snakemake.params.res
res_ratio = snakemake.params.res_ratio
width = snakemake.params.width
cpus = snakemake.threads

outfile = str(snakemake.output.matrices)
outfile_ratio = str(snakemake.output.ratio)
outfile_ratio_rif = str(snakemake.output.ratio_rif)
# out_zoom = str(snakemake.output.zoom)
out_zoom_ratio = str(snakemake.output.zoom_ratio)

# Make sure output directory exists.
os.makedirs(str((snakemake.params.outdir)), exist_ok=True)

# Import matrices
subsample = 11675151
files = [mat_t7_ara_file, mat_t7_novo_file, mat_t7_ara_rif_file, mat_t7_novo_rif_file]
hic = [0, 0, 0, 0]
for i, file in enumerate(files):
    clr = cooler.Cooler(f"{file}::/resolutions/{res}")
    mat = clr.matrix(balance=False, sparse=True)[:, :]
    mat = hcs.subsample_contacts(mat, subsample)
    mat = hcs.normalize_sparse(mat, norm='ICE', iterations=100, n_mad=10)
    hic[i] = mat.toarray()

# Ratio
subsample = 10251593
hic_ratio = [0, 0, 0, 0]
for i, file in enumerate(files):
    clr = cooler.Cooler(f"{file}::/resolutions/{res_ratio}")
    mat = clr.matrix(balance=False, sparse=True)[:, :]
    print(mat.sum())
    mat = hcs.subsample_contacts(mat, subsample)
    mat = hcs.normalize_sparse(mat, norm='ICE', iterations=100, n_mad=10)
    hic_ratio[i] = mat.toarray()

ratio = np.log10(hic_ratio[1]) - np.log10(hic_ratio[0])
ratio[np.isnan(ratio)] =  0
ratio[ratio == np.inf] =  0
ratio[ratio == -np.inf] =  0
ratio_rif = np.log10(hic_ratio[3]) - np.log10(hic_ratio[2])
ratio_rif[np.isnan(ratio_rif)] =  0
ratio_rif[ratio_rif == np.inf] =  0
ratio_rif[ratio_rif == -np.inf] =  0

###
# PLOT MATRICES
###

title = ["novo- Rif-", "novo+ Rif-", "novo- Rif+", "novo+ Rif+"]
ratios = [ratio, ratio_rif]
start = 0
end = len(hic[0])
end_ratio = end * res // res_ratio 

fig, ax = plt.subplots(3, 2, figsize=(10, 10))
for i in range(2):
    for j in range(2):
        ax[i, j].imshow(
            hic[i + 2 * j][start:end, start:end],
            vmin=0,
            vmax=0.0012,
            cmap="Reds",
            extent=(start, end, end, start),
        )
    ax[i, j].set_title(title[i + 2 * j])
    ax[2, i].imshow(
        ratios[i][start:end_ratio, start:end_ratio],
        vmin=-2,
        vmax=2,
        cmap="seismic",
        extent=(start, end, end, start),
    )
    ax[2, i].set_title("novo+/novo-")

plt.savefig(outfile)
plt.close()

###
# Plot only ratio.
###

fig, ax = plt.subplots(1, 1, figsize=(10, 10))
ax.imshow(
    ratios[0][start:end_ratio, start:end_ratio],
    vmin=-2,
    vmax=2,
    cmap="seismic",
    extent=(start, end, end, start),
)
ax.set_title("novo+/novo-")

plt.savefig(outfile_ratio, dpi=300)
plt.close()

fig, ax = plt.subplots(1, 1, figsize=(10, 10))
ax.imshow(
    ratios[1][start:end_ratio, start:end_ratio],
    vmin=-2,
    vmax=2,
    cmap="seismic",
    extent=(start, end, end, start),
)
ax.set_title("novo+ Rif+/novo- Rif+")

plt.savefig(outfile_ratio_rif, dpi=300)
plt.close()

###
# Plot zoom without RNAseq.
###

# Rotate matrices.
# hic_rot = copy.copy(hic)
# for i in range(len(hic)):
#     hic_rot[i] = bch.interpolate_white_lines(hic_rot[i])
#     hic_rot[i][np.isnan(hic_rot[i])] = 0
#     hic_rot[i] = ndimage.rotate(hic_rot[i], 45, reshape=True)

# fig, ax = plt.subplots(
#     1, 4, figsize=(20, 10)
# )

# start = 300_000
# end = 480_000

# row1 = len(hic_rot[0]) // 2 - int(np.sqrt(2) * (width / res))
# row2 = len(hic_rot[0]) // 2 + int(np.sqrt(2) * (width / res))
# col1 = int((start // res) * np.sqrt(2))
# col2 = int((end // res) * np.sqrt(2))

# for i in range(len(hic)):
#     # Hicmap
#     im = ax[i].imshow(
#         hic_rot[i][row1:row2, col1:col2],
#         cmap="Reds",
#         vmin=0,
#         vmax=0.002,
#         extent=(col1 / np.sqrt(2), col2 / np.sqrt(2), width // res, -width // res),
#     )
#     ax[i].get_xaxis().set_visible(False)
#     ax[i].set_title(title[i], size=18)
#     if i == 0:
#         ax[i].tick_params(axis="both", labelsize=16)
#         ax[i].set_ylabel("Genomic coordinates (kb)", size=16)
#     else:
#         ax[i].get_yaxis().set_visible(False)

# # Colorbar
# cbar = fig.colorbar(
#     im,
#     ax=ax.ravel().tolist(),
#     shrink=0.5,
#     anchor=(1., 0.5)
# )
# cbar.ax.tick_params(labelsize=16)

# # Adjust space between plot
# plt.subplots_adjust(wspace=0.4, hspace=0.0, left=0.05, right=0.88)
# plt.savefig(out_zoom)
# plt.close()

###
# Plot zoom ratio.
###

start, end = 0, 750
end_ratio = end * res // res_ratio 

fig, ax = plt.subplots(3, 2, figsize=(10, 10))
for i in range(2):
    for j in range(2):
        ax[i, j].imshow(
            hic[i + 2 * j][start:end, start:end],
            vmin=0,
            vmax=0.0012,
            cmap="Reds",
            extent=(start, end, end, start),
        )
    ax[i, j].set_title(title[i + 2 * j])
for i in range(2):
    ax[2, i].imshow(
        ratios[i][start:end_ratio, start:end_ratio],
        vmin=-1,
        vmax=1,
        cmap="seismic",
        extent=(start, end, end, start),
    )
    ax[2, i].set_title("novo+/novo-")

plt.savefig(out_zoom_ratio)
plt.close() 
