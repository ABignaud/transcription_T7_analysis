#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import bacchus.hic as bch
import bacchus.io as bcio
import cooler
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.stats as st

# Import snakemake.
mat_t7 = snakemake.input.mat_t7
mat_t7_2p100 = snakemake.input.mat_t7_2p100
mat_t7_2p60 = snakemake.input.mat_t7_2p60
mat_t7_2p100_dv = snakemake.input.mat_t7_2p100_dv
mat_t7_2p60_dv = snakemake.input.mat_t7_2p60_dv
mat_t7_2p100_cv = snakemake.input.mat_t7_2p100_cv
rna_t7 = snakemake.input.rna_t7
rna_t7_2p100 = snakemake.input.rna_t7_2p100
rna_t7_2p60 = snakemake.input.rna_t7_2p60
rna_t7_2p100_dv = snakemake.input.rna_t7_2p100_dv
rna_t7_2p60_dv = snakemake.input.rna_t7_2p60_dv
rna_t7_2p100_cv = snakemake.input.rna_t7_2p100_cv
chip_RNA_t7 = snakemake.input.chip_RNA_t7
chip_RNA_t7_2p100 = snakemake.input.chip_RNA_t7_2p100
chip_RNA_t7_2p60 = snakemake.input.chip_RNA_t7_2p60
chip_RNA_t7_2p100_dv = snakemake.input.chip_RNA_t7_2p100_dv
chip_RNA_t7_2p60_dv = snakemake.input.chip_RNA_t7_2p60_dv
chip_RNA_t7_2p100_cv = snakemake.input.chip_RNA_t7_2p100_cv
chip_gapr_t7 = snakemake.input.chip_gapr_t7
chip_gapr_t7_2p100 = snakemake.input.chip_gapr_t7_2p100
chip_gapr_t7_2p60 = snakemake.input.chip_gapr_t7_2p60
chip_gapr_t7_2p100_dv = snakemake.input.chip_gapr_t7_2p100_dv
chip_gapr_t7_2p60_dv = snakemake.input.chip_gapr_t7_2p60_dv
chip_gapr_t7_2p100_cv = snakemake.input.chip_gapr_t7_2p100_cv
chip_gapr_control = snakemake.input.chip_gapr_control
cmap = snakemake.params.cmap
res = snakemake.params.res
tracks_res = snakemake.params.tracks_res 
outdir = snakemake.params.outdir 
out_plot = str(snakemake.output.plot)
out_txt = str(snakemake.output.txt)

# Make sure output directory exists.
os.makedirs(outdir, exist_ok=True)

# Import matrices
mat_t7 = cooler.Cooler(f"{mat_t7}::/resolutions/{res}").matrix(
    balance=True, sparse=False
)[:]
mat_t7[np.isnan(mat_t7)] = 0
mat_t7_2p100 = cooler.Cooler(f"{mat_t7_2p100}::/resolutions/{res}").matrix(
    balance=True, sparse=False
)[:]
mat_t7_2p100[np.isnan(mat_t7_2p100)] = 0
mat_t7_2p60 = cooler.Cooler(f"{mat_t7_2p60}::/resolutions/{res}").matrix(
    balance=True, sparse=False
)[:]
mat_t7_2p60[np.isnan(mat_t7_2p60)] = 0
mat_t7_2p100_dv = cooler.Cooler(f"{mat_t7_2p100_dv}::/resolutions/{res}").matrix(
    balance=True, sparse=False
)[:]
mat_t7_2p100_dv[np.isnan(mat_t7_2p100_dv)] = 0
mat_t7_2p60_dv = cooler.Cooler(f"{mat_t7_2p60_dv}::/resolutions/{res}").matrix(
    balance=True, sparse=False
)[:]
mat_t7_2p60_dv[np.isnan(mat_t7_2p60_dv)] = 0
mat_t7_2p100_cv = cooler.Cooler(f"{mat_t7_2p100_cv}::/resolutions/{res}").matrix(
    balance=True, sparse=False
)[:]
mat_t7_2p100_cv[np.isnan(mat_t7_2p100_cv)] = 0

# Import RNA tracks
rna_t7, _ = bcio.extract_big_wig(rna_t7, tracks_res)
rna_t7_2p100, _ = bcio.extract_big_wig(rna_t7_2p100, tracks_res)
rna_t7_2p60, _ = bcio.extract_big_wig(rna_t7_2p60, tracks_res)
rna_t7_2p100_dv, _ = bcio.extract_big_wig(rna_t7_2p100_dv, tracks_res)
rna_t7_2p60_dv, _ = bcio.extract_big_wig(rna_t7_2p60_dv, tracks_res)
rna_t7_2p100_cv, _ = bcio.extract_big_wig(rna_t7_2p100_cv, tracks_res)

# Import chip RNA pol T7
chip_RNA_t7, _ = bcio.extract_big_wig(chip_RNA_t7, res)
chip_RNA_t7_2p100, _ = bcio.extract_big_wig(chip_RNA_t7_2p100, res)
chip_RNA_t7_2p60, _ = bcio.extract_big_wig(chip_RNA_t7_2p60, res)
chip_RNA_t7_2p100_dv, _ = bcio.extract_big_wig(chip_RNA_t7_2p100_dv, res)
chip_RNA_t7_2p60_dv, _ = bcio.extract_big_wig(chip_RNA_t7_2p60_dv, res)
chip_RNA_t7_2p100_cv, _ = bcio.extract_big_wig(chip_RNA_t7_2p100_cv, res)

# Import chip GapR
chip_gapr_t7, _ = bcio.extract_big_wig(chip_gapr_t7, res)
chip_gapr_t7_2p100, _ = bcio.extract_big_wig(chip_gapr_t7_2p100, res)
chip_gapr_t7_2p60, _ = bcio.extract_big_wig(chip_gapr_t7_2p60, res)
chip_gapr_t7_2p100_dv, _ = bcio.extract_big_wig(chip_gapr_t7_2p100_dv, res)
chip_gapr_t7_2p60_dv, _ = bcio.extract_big_wig(chip_gapr_t7_2p60_dv, res)
chip_gapr_t7_2p100_cv, _ = bcio.extract_big_wig(chip_gapr_t7_2p100_cv, res)

# Control GapR
chip_gapr_control, _ = bcio.extract_big_wig(chip_gapr_control, res)

# Compute the log fold change between signal
chip_gapr_t7 = np.log2(chip_gapr_t7) - np.log2(chip_gapr_control) 
chip_gapr_t7_2p100 = np.log2(chip_gapr_t7_2p100) - np.log2(chip_gapr_control)
chip_gapr_t7_2p60 = np.log2(chip_gapr_t7_2p60) - np.log2(chip_gapr_control)
chip_gapr_t7_2p100_dv = np.log2(chip_gapr_t7_2p100_dv) - np.log2(chip_gapr_control)
chip_gapr_t7_2p60_dv = np.log2(chip_gapr_t7_2p60_dv) - np.log2(chip_gapr_control)
chip_gapr_t7_2p100_cv = np.log2(chip_gapr_t7_2p100_cv) - np.log2(chip_gapr_control) 

# Compute HiC signal
hic_t7 = bch.compute_hic_signal(mat_t7, 1000, 0, 5000)
hic_t7_2p100 = bch.compute_hic_signal(mat_t7_2p100, 1000, 0, 5000)
hic_t7_2p60 = bch.compute_hic_signal(mat_t7_2p60, 1000, 0, 5000)
hic_t7_2p100_dv = bch.compute_hic_signal(mat_t7_2p100_dv, 1000, 0, 5000)
hic_t7_2p60_dv = bch.compute_hic_signal(mat_t7_2p60_dv, 1000, 0, 5000)
hic_t7_2p100_cv = bch.compute_hic_signal(mat_t7_2p100_cv, 1000, 0, 5000)

list_rna = [
    rna_t7,
    rna_t7_2p100,
    rna_t7_2p60,
    rna_t7_2p100_dv,
    rna_t7_2p60_dv,
    rna_t7_2p100_cv,
]
list_chip_T7 = [
    chip_RNA_t7,
    chip_RNA_t7_2p100,
    chip_RNA_t7_2p60,
    chip_RNA_t7_2p100_dv,
    chip_RNA_t7_2p60_dv,
    chip_RNA_t7_2p100_cv,
]
list_chip_gapr = [
    chip_gapr_t7,
    chip_gapr_t7_2p100,
    chip_gapr_t7_2p60,
    chip_gapr_t7_2p100_dv,
    chip_gapr_t7_2p60_dv,
    chip_gapr_t7_2p100_cv,
]
list_hic_signal = [
    hic_t7,
    hic_t7_2p100,
    hic_t7_2p60,
    hic_t7_2p100_dv,
    hic_t7_2p60_dv,
    hic_t7_2p100_cv,
]
list_hic = [
    mat_t7,
    mat_t7_2p100,
    mat_t7_2p60,
    mat_t7_2p100_dv,
    mat_t7_2p60_dv,
    mat_t7_2p100_cv,
]

def z_transform(values):
    return (values - np.nanmean(values) / np.nanstd(values))

# for i in range(6):
#     list_rna[i] = z_transform(list_rna[i])
#     list_chip_T7[i] = z_transform(list_chip_T7[i])
#     list_chip_gapr[i] = z_transform(list_chip_gapr[i])
#     list_hic_signal[i] = z_transform(list_hic_signal[i])

# Make the plots 
start, end = 0, 750
binning = 10

fig, ax = plt.subplots(4,6 , figsize=(20,15), sharex=True, gridspec_kw={'height_ratios': [8, 4, 4, 4]})

# Define max values for x_lim
max_rna = 15.5
min_gapr = -4

# Define color
col = ["k", "#fdbf6f", "#1f78b4", '#e31a1c']

for i in range(6):
    # Hicmap
    im = ax[0, i].imshow(
        list_hic[i][start:end, start:end],
        cmap="Reds",
        vmax=np.nanpercentile(list_hic[0], 99.9),
        extent=(start, end, end, start),
    )
    ax[0, i].get_xaxis().set_visible(False)
    ax[0, i].tick_params(axis='both', labelsize=16)
    # ax[0, i].set_ylabel("Genomic coordinates (kb)", size=16)
    
    # Settings plot 1, 2
    # ax[1, i].set_xlabel("Genomic coordinates (kb)", size=16)
    ax[1, 0].set_ylabel("Transcription (CPM)", size=16)
    ax[1, i].spines['top'].set_visible(False)
    ax[1, i].spines['right'].set_visible(False)
    ax[1, i].set_ylim(-1, 3000)
    ax[1, i].tick_params(axis='both', labelsize=16, color=col[0], labelcolor=col[0])
    # ax[2, i].set_xlabel("Genomic coordinates (kb)", size=16)
    ax[2, 0].set_ylabel("T7 RNA pol signal (CPM)", size=16)
    ax[2, i].spines['top'].set_visible(False)
    ax[2, i].spines['right'].set_visible(False)
    ax[2, i].set_ylim(-0.05, 1)
    ax[2, i].tick_params(axis='both', labelsize=16, color=col[0], labelcolor=col[0])
    ax[3, i].set_xlabel("Genomic coordinates (kb)", size=16)
    ax[3, 0].set_ylabel("GapR signal (log2 fold change)", size=16)
    ax[3, i].spines['top'].set_visible(False)
    ax[3, i].spines['right'].set_visible(False)
    ax[3, i].set_ylim(-0.05, 1)
    ax[3, i].tick_params(axis='both', labelsize=16, color=col[0], labelcolor=col[0])
    
    # RNAseq
    if list_rna[i].all() != None:
        p1 = ax[1, i].fill_between(
            np.arange(start, end, 1 / binning),
            list_rna[i][int(start * binning):int(end * binning)],
            color=col[0],
            label="RNAseq",
            linewidth=1,
        )
    
    # HiC signal
    p2 = ax[2, i].plot(
        np.arange(start, end, 1),
        list_hic_signal[i][start:end],
        color=col[3],
        label="HiC signal",
        linewidth=1,
        alpha=.7,
    )
    p2 = ax[3, i].plot(
        np.arange(start, end, 1),
        list_hic_signal[i][start:end],
        color=col[3],
        label="HiC signal",
        linewidth=1,
        alpha=.7,
    )
    
    # T7 ChIPseq
    ax3 = ax[2, i].twinx()
    if list_chip_T7[i].all() != None:
        p3 = ax3.plot(
            np.arange(start, end, 1),
            list_chip_T7[i][start:end],
            color=col[2],
            label="T7 RNA pol ChIPseq",
            alpha=.7,
            linewidth=1,
        )
        ax3.set_ylim(-100, 6000)
    
    # GapR ChIPseq
    ax4 = ax[3, i].twinx()
    if list_chip_gapr[i].all() != None:
        p4 = ax4.plot(
            np.arange(start, end, 1),
            list_chip_gapr[i][start:end],
            color=col[1],
            label="GapR ChIPseq",
            alpha=1,
            linewidth=1,
        )
        # ax[3, i].set_ylim(-7.5, 10)
        ax4.set_ylim(-1.5, 1.5)
    
# Colorbar
cbar = fig.colorbar(im, ax=ax.ravel().tolist(), shrink=.2, anchor=(1.2, .89))
cbar.ax.tick_params(labelsize=16)

# Legend
leg = fig.legend(
    labels=["RNAseq", "GapR ChIPseq", "T7 RNA pol ChIPseq", "HiC signal"],
    loc=(0.9, 0.25),
)
for i in range(4):
    leg.legendHandles[i].set_color(col[i])

# Adjust space between plot
plt.subplots_adjust(wspace=.4, hspace=.2, left=0.05, right=0.88)
plt.savefig(out_plot, dpi=250)


# Correlation between tracks
start=220 #kb
end=580 #kb

def get_binning(values, size=10):
    n = len(values) // size
    new_values = np.zeros(n + 1)
    for i in range(n):
        new_values[i] = np.nanmean(values[i * size:(i + 1) * size])
    new_values[n] = np.nanmean(values[n * size:])
    new_values
    return new_values

for i in range(6):
    list_rna[i] = get_binning(list_rna[i])
    # list_chip_T7[i] = get_binning(list_chip_T7[i])
    # list_chip_gapr[i] = get_binning(list_chip_gapr[i])


with open(out_txt, 'w') as out:
    out.write('Correlation between signal:\n')
    label=["T7", "T7_CL_100", "T7_CL_60", "T7_DV_100", "T7_DV_60", "T7_CV"]
    for i in range(6):
    #     out.write(f"{label[i]};RNA/T7 Pearson correlation: {st.pearsonr(list_rna[i][start:end], list_chip_T7[i][start:end])[0]:.2f}\n")
    #     out.write(f"{label[i]};RNA/GapR Pearson correlation: {st.pearsonr(list_rna[i][start:end], list_chip_gapr[i][start:end])[0]:.2f}\n")
    #     out.write(f"{label[i]};RNA/HiC Pearson correlation: {st.pearsonr(list_rna[i][start:end], list_hic_signal[i][start:end])[0]:.2f}\n")
    #     out.write(f"{label[i]};T7/GapR Pearson correlation: {st.pearsonr(list_chip_T7[i][start:end], list_chip_gapr[i][start:end])[0]:.2f}\n")
    #     out.write(f"{label[i]};T7/HiC Pearson correlation: {st.pearsonr(list_chip_T7[i][start:end], list_hic_signal[i][start:end])[0]:.2f}\n")
    #     out.write(f"{label[i]};GapR/HiC Pearson correlation: {st.pearsonr(list_chip_gapr[i][start:end], list_hic_signal[i][start:end])[0]:.2f}\n")
        out.write(f"{label[i]};RNA/T7 Spearmann correlation: {st.spearmanr(list_rna[i][start:end], list_chip_T7[i][start:end])[0]:.2f}\n")    
        out.write(f"{label[i]};RNA/GapR Spearmann correlation: {st.spearmanr(list_rna[i][start:end], list_chip_gapr[i][start:end])[0]:.2f}\n")
        out.write(f"{label[i]};RNA/HiC Spearmann correlation: {st.spearmanr(list_rna[i][start:end], list_hic_signal[i][start:end])[0]:.2f}\n")
        out.write(f"{label[i]};T7/GapR Spearmann correlation: {st.spearmanr(list_chip_T7[i][start:end], list_chip_gapr[i][start:end])[0]:.2f}\n")
        out.write(f"{label[i]};T7/HiC Spearmann correlation: {st.spearmanr(list_chip_T7[i][start:end], list_hic_signal[i][start:end])[0]:.2f}\n")  
        out.write(f"{label[i]};GapR/HiC Spearmann correlation: {st.spearmanr(list_chip_gapr[i][start:end], list_hic_signal[i][start:end])[0]:.2f}\n")
