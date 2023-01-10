#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import bacchus.directional as bcd
import cooler
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.sparse as sp

mat_file = snakemake.input.mat
res = snakemake.params.res
cmap = snakemake.params.cmap
macro_output = str(snakemake.output.macrodomain)
CID_output = str(snakemake.output.cid)
out_file = str(snakemake.output.compare)

# Create outdir if necessary.
os.makedirs(str(snakemake.params.outdir), exist_ok=True)

mat = cooler.Cooler(f"{mat_file}::/resolutions/{res}").matrix(
    balance=True, sparse=True
)[:]

# Positions of the CIDs borders from Lioy et al., Cell, 2018.
lioy = [
    40,
    145,
    220,
    440,
    510,
    760,
    980,
    1145,
    1195,
    1300,
    1585,
    1670,
    1795,
    1845,
    2100,
    2255,
    2380,
    2525,
    2720,
    2895,
    2990,
    3165,
    3295,
    3430,
    3645,
    3800,
    3935,
    4030,
    4160,
    4205,
    4465,
]

# Length of E coli chromosome in kb
n = 4642

# Comput min max size
max_cid = n - lioy[-1] + lioy[0]
min_cid = n - lioy[-1] + lioy[0]
prev_cid = lioy[0]
for i in lioy[1:]:
    size = i - prev_cid
    prev_cid = i
    if size > max_cid:
        max_cid = size
    if size < min_cid:
        min_cid = size

print(f"Longest CID- Lioy: {max_cid}kb")
print(f"Shortest CID - Lioy: {min_cid}kb")

# Compute macrodomains on our matrix based on directional index.
di_macro = bcd.directional_index(mat, 80)
plt.subplots(figsize=(12, 1))
plt.fill_between(
    x=np.arange(0, len(di_macro) * 5, 5),
    y1=0,
    y2=di_macro,
    where=di_macro > 0,
    color="g",
)
plt.fill_between(
    x=np.arange(0, len(di_macro) * 5, 5),
    y1=0,
    y2=di_macro,
    where=di_macro <= 0,
    color="r",
)
plt.ylim(-2, 2)
plt.xlim(0, n)
plt.savefig(macro_output)

borders_macro = bcd.di_borders(di_macro)
print("Numbers of CIDs:", len(borders_macro))

# Compute CIDs on our matrix based on directional index.
di_CIDs = bcd.directional_index(mat, 20)
plt.subplots(figsize=(12, 2))
plt.fill_between(
    x=np.arange(0, len(di_CIDs) * 5, 5),
    y1=0,
    y2=di_CIDs,
    where=di_CIDs > 0,
    color="g",
)
plt.fill_between(
    x=np.arange(0, len(di_CIDs) * 5, 5),
    y1=0,
    y2=di_CIDs,
    where=di_CIDs <= 0,
    color="r",
)
plt.ylim(-2, 2)
plt.xlim(0, n)
plt.savefig(CID_output)

borders_CIDs = bcd.di_borders(di_CIDs)
print("Numbers of CIDs:", len(borders_CIDs))

from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# Plot the CID on a matrix at 5kb with their rspectives positions marked with
# stars.
axis = "kb"
start = 0
title = "WT E. coli domains - binning 5kb"
dpi = 300
vmax = 99

fig, ax = plt.subplots(
    3,
    1,
    figsize=(11, 13),
    dpi=dpi,
    gridspec_kw={"height_ratios": [10, 1, 1]},
    sharex=True,
)

# Axis values
scaling_factor = res // 1000
end = n // scaling_factor

# Display plots
im = ax[0].imshow(
    mat.toarray()[start:end, start:end],
    cmap=cmap,
    vmin=0,
    vmax=np.percentile(mat.toarray(), vmax),
    extent=(
        start * scaling_factor,
        end * scaling_factor,
        end * scaling_factor,
        start * scaling_factor,
    ),
)

# Legend
ax[2].set_xlabel(f"Genomic coordinates ({axis:s})", fontsize=16)
ax[0].set_ylabel(f"Genomic coordinates ({axis:s})", fontsize=16)
ax[0].tick_params(axis="both", labelsize=16)
ax[1].tick_params(axis="both", labelsize=16)
ax[2].tick_params(axis="both", labelsize=16)
ax[2].tick_params(axis="x", which="major", pad=15)
ax[1].set_ylabel("DI\n(400kb)", fontsize=16)
ax[2].set_ylabel("DI\n(100kb)", fontsize=16)
ax[0].set_title(title, size=18)
ax[1].set_title("Macrodomains", size=16)
ax[2].set_title("CIDs", size=16)

# Colorbar
cbar = plt.colorbar(im, ax=ax.ravel().tolist(), shrink=0.33, anchor=(0, 0.7))
cbar.ax.tick_params(labelsize=16)

# Add DI plots
ax[1].fill_between(
    x=np.arange(0, len(di_macro) * 5, 5),
    y1=0,
    y2=di_macro,
    where=di_macro > 0,
    color="#33a02c",
)
ax[1].fill_between(
    x=np.arange(0, len(di_macro) * 5, 5),
    y1=0,
    y2=di_macro,
    where=di_macro <= 0,
    color="#e31a1c",
)
ax[1].set_ylim(-2, 2)
ax[1].set_xlim(start, end * scaling_factor)
for i in borders_macro:
    if i > start and i < end:
        ax[1].text(x=(i * 5) - 25, y=-3, s="*", color="k", fontweight="bold")

ax[2].fill_between(
    x=np.arange(0, len(di_CIDs) * 5, 5),
    y1=0,
    y2=di_CIDs,
    where=di_CIDs >= 0,
    color="#33a02c",
    interpolate=True,
)
ax[2].fill_between(
    x=np.arange(0, len(di_CIDs) * 5, 5),
    y1=0,
    y2=di_CIDs,
    where=di_CIDs <= 0,
    color="#e31a1c",
    interpolate=True,
)
ax[2].set_ylim(-2, 2)
ax[2].set_xlim(start * scaling_factor, end * scaling_factor)
for i in borders_CIDs:
    if i > start and i < end:
        ax[0].axvline(x=i * 5, linestyle="dashed", color="k")
        ax[2].text(x=(i * 5) - 25, y=-3, s="*", color="k", fontweight="bold")
for i in [133, 176, 394, 713, 872]:
    if i > start and i < end:
        ax[2].text(
            x=(i * 5) - 25, y=-3, s="*", color="#33a02c", fontweight="bold"
        )
for i in [102, 239, 317, 369, 420, 505, 633, 729, 832]:
    if i > start and i < end:
        ax[2].text(
            x=(i * 5) - 25, y=-3, s="*", color="#e31a1c", fontweight="bold"
        )
# Add legend of stars
legend_elements = [
    Line2D(
        [],
        [],
        color="k",
        marker="*",
        linestyle="None",
        markersize=10,
        label="Conserved borders",
    ),
    Line2D(
        [],
        [],
        color="#33a02c",
        marker="*",
        linestyle="None",
        markersize=10,
        label="New borders",
    ),
    Line2D(
        [],
        [],
        color="#e31a1c",
        marker="*",
        linestyle="None",
        markersize=10,
        label="Missing borders",
    ),
]
ax[2].legend(handles=legend_elements, bbox_to_anchor=(1.3, 1.05))

# Savefig
if out_file is not None:
    plt.savefig(out_file, dpi=dpi)
