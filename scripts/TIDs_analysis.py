#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import bacchus.blob as bcb
import bacchus.hic as bch
import bacchus.io as bcio
import numpy as np
import copy
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import scipy
import scipy.stats as st
import seaborn as sns
import cooler

mat_file = snakemake.input.mat
rna_file = snakemake.input.rna
cov_file = snakemake.input.cov
gc_file = snakemake.input.gc
frags_file = snakemake.input.frags
epod_file = snakemake.input.EPODs
out_bed = str(snakemake.output.bed)
text_file = str(snakemake.output.text_file)
gc_plot = str(snakemake.output.gc)
cov_plot = str(snakemake.output.cov)
RS_plot = str(snakemake.output.RS)
rna_plot = str(snakemake.output.rna)
dist_bar = str(snakemake.output.dist_bar)
dist_line = str(snakemake.output.dist_line)
venn_plot = str(snakemake.output.venn_plot)

# Create outdir if necessary.
os.makedirs(str(snakemake.params.outdir), exist_ok=True)


def write_stats(text_file, text, writing_type="a"):
    with open(text_file, writing_type) as out:
        out.write(f"{text}\n")


# Import data
mat = cooler.Cooler(f"{mat_file}::/resolutions/500").matrix(balance=True)[:]
mat[np.isnan(mat)] = 0
rna, _ = bcio.extract_big_wig(rna_file, ztransform=False)
rna500, _ = bcio.extract_big_wig(rna_file, binning=500, ztransform=False)
rna500 = np.log10(rna500[:-2] + 1e-10)
cov, _ = bcio.extract_big_wig(cov_file, binning=500)
cov = cov[:-2]
gc_content = pd.read_csv(gc_file, sep="\t", header=None).iloc[:, 1]
frags = pd.read_csv(frags_file, sep="\t")
epod_data = pd.read_csv(epod_file, sep="\t", header=None)
mask = bch.mask_white_line(mat)
mat[mask] = np.nan
mat[:, mask] = np.nan

# Save TIDs positions.
blobs, blob_score = bcb.find_blobs(
    mat, size=5, n_mads=10, refine=1 / 3, rna=rna
)
with open(out_bed, "w") as out:
    for blob in blobs:
        out.write(f"blob\t{blob.start * 500}\t{blob.end * 500}\n")

# Save stats on TIDs.
write_stats(
    text_file, f"Numbers of blobs without interpolation: {len(blobs)}", "w"
)
write_stats(
    text_file,
    f"Cumulative size of the blobs without interpolation: {np.sum([x.size for x in blobs]) * 0.5}kb",
)

# Mask bins within a TIDs or not.
blob_mask = np.repeat(False, len(mat))
for i in blobs:
    for j in range(i.start, i.end):
        blob_mask[j] = True
blob_mask_r = np.repeat(False, len(mat))
for i, v in enumerate(blob_mask):
    if v:
        blob_mask_r[i] = False
    else:
        blob_mask_r[i] = True

# GC content
gc_content = np.array(gc_content[np.arange(0, len(gc_content), 5)])
gc_blob = gc_content[blob_mask]
gc_other = gc_content[blob_mask_r]
data = {"GC": gc_content, "blobs": blob_mask}
data = pd.DataFrame(data)
sns.violinplot(x="blobs", y="GC", data=data, palette="tab10")
plt.xlabel("Blobs", size=14)
plt.ylabel("GC content %", size=14)
plt.title("GC content distribution in blobs", size=16)
plt.xticks(size=12)
plt.yticks(size=12)
plt.savefig(gc_plot, bbox_inches="tight")
plt.close()

write_stats(text_file, f"Blob GC content: {np.mean(gc_blob)}")
write_stats(text_file, f"Other GC content: {np.mean(gc_other)}")
write_stats(text_file, f"Global GC content: {np.mean(gc_content)}")
write_stats(text_file, f"T-test blob-other: {st.ttest_ind(gc_blob, gc_other)}")
write_stats(text_file, f"T-test blob-all: {st.ttest_ind(gc_blob, gc_content)}")

# Restriction sites
start_pos = frags.start_pos[1:]
frags = np.zeros(len(gc_content))
for i in start_pos:
    frags[i // 500] += 1
data["frags"] = frags

sns.boxplot(x="blobs", y="frags", data=data, palette="tab10")
plt.xlabel("Blobs", size=14)
plt.ylabel("Restriction site", size=14)
plt.title("Restriction sites distribution in blobs", size=16)
plt.xticks(size=12)
plt.yticks(size=12)
plt.savefig(RS_plot, bbox_inches="tight")
plt.close()

write_stats(text_file, f"Blob RS: {np.mean(frags[blob_mask])}")
write_stats(text_file, f"Other RS: {np.mean(frags[blob_mask_r])}")
write_stats(text_file, f"Global RS: {np.mean(frags)}")
write_stats(
    text_file,
    f"T-test blob-other: {st.ttest_ind(frags[blob_mask], frags[blob_mask_r])}",
)
write_stats(
    text_file, f"T-test blob-all: {st.ttest_ind(frags[blob_mask], frags)}"
)

# Coverage
data["cov"] = cov

sns.violinplot(x="blobs", y="cov", data=data, palette="tab10")
plt.xlabel("Blobs", size=14)
plt.ylabel("HiC Coverage (CPM)", size=14)
plt.title("Coverage distribution in blobs", size=16)
plt.xticks(size=12)
plt.yticks(size=12)
plt.savefig(cov_plot, bbox_inches="tight")
plt.close()

write_stats(text_file, f"Blob coverage: {np.mean(cov[blob_mask])}")
write_stats(text_file, f"Other coverage: {np.mean(cov[blob_mask_r])}")
write_stats(text_file, f"Global coverage: {np.mean(cov)}")
write_stats(
    text_file,
    f"T-test blob-other: {st.ttest_ind(cov[blob_mask], cov[blob_mask_r])}",
)
write_stats(text_file, f"T-test blob-all: {st.ttest_ind(cov[blob_mask], cov)}")

# Transcription
data["rna"] = rna500

sns.violinplot(x="blobs", y="rna", data=data, palette="tab10")
plt.xlabel("Blobs", size=14)
plt.ylabel("Transcription (CPM)", size=14)
plt.title("Transcription distribution in blobs", size=16)
plt.xticks(size=12)
plt.yticks(size=12)
plt.ylim(-8, 6)
plt.savefig(rna_plot, bbox_inches="tight")
plt.close()

write_stats(text_file, f"Blob transcription: {np.mean(rna500[blob_mask])}")
write_stats(text_file, f"Other transcription: {np.mean(rna500[blob_mask_r])}")
write_stats(text_file, f"Global transcription: {np.mean(rna500)}")
write_stats(
    text_file,
    f"T-test blob-other: {st.ttest_ind(rna500[blob_mask], rna500[blob_mask_r])}",
)
write_stats(
    text_file, f"T-test blob-all: {st.ttest_ind(rna500[blob_mask], rna500)}"
)

# TIDs distribution

tid_data = pd.read_csv(out_bed, sep="\t", header=None)
ori = (3925744 + 3925975) // 2
ter = 1590754
macrodomain = [x * 5000 for x in [45, 231, 393, 546, 685, 842]]
x = np.arange(0, 4_641_652, 500)
y = np.zeros(len(x))
for i in tid_data.index:
    start = tid_data.loc[i, 1]
    end = tid_data.loc[i, 2]
    for k in np.arange(start, end, 500):
        y[k // 500] = 1

# Bar plot
fig, ax = plt.subplots(figsize=(20, 5))
ax.fill_between(x / 1_000_000, y, color="k", alpha=1)
plt.xticks(size=16)
plt.xlabel("Genomics coordinates (Mb)", size=18)
ax.yaxis.set_visible(False)
ax.set_title("TIDs distribution", size=20)
for i in macrodomain[:-1]:
    plt.axvline(
        x=i / 1_000_000, c="green", linestyle="dashed", alpha=0.5, linewidth=3
    )
plt.axvline(
    x=macrodomain[-1] / 1_000_000,
    c="green",
    linestyle="dashed",
    alpha=0.5,
    label="macrodomain",
    linewidth=3,
)
plt.axvline(
    x=ori / 1_000_000,
    c="r",
    linestyle="dashed",
    alpha=0.5,
    label="ori",
    linewidth=3,
)
plt.axvline(
    x=ter / 1_000_000,
    c="blue",
    linestyle="dashed",
    alpha=0.5,
    label="ter",
    linewidth=3,
)
plt.legend()
plt.savefig(dist_bar)
plt.close()

# Line plot
bins = np.arange(0, 4_641_652, 50_000)
hist, x_bins = np.histogram(x, bins, weights=y)
fig, ax = plt.subplots(figsize=(15, 5))
ax.plot(x_bins[:-1] / 1_000_000, hist * 1, color="k")
ax.set_xlabel("Genomic coordinates (Mb)", size=16)
ax.set_ylabel("Ratio of TID in a 10kb bin (%)", size=16)
ax.tick_params(labelsize=16)
ax.set_title("TIDs distribution", size=20)
ax.set_xlim(0, 4_641_652 / 1_000_000)
for i in macrodomain[:-1]:
    plt.axvline(x=i / 1_000_000, c="k", linestyle="dashed", alpha=0.5)
plt.axvline(
    x=macrodomain[-1] / 1_000_000,
    c="k",
    linestyle="dashed",
    alpha=0.5,
    label="macrodomain",
)
plt.axvline(
    x=ori / 1_000_000, c="r", linestyle="dashed", alpha=0.5, label="ori"
)
plt.axvline(
    x=ter / 1_000_000, c="blue", linestyle="dashed", alpha=0.5, label="ter"
)
plt.legend()
plt.savefig(dist_line)
plt.close()

# Venn Diagramm

x = np.arange(0, 4_641_652, 5)
y_tid = np.zeros(len(x))
y_epod = np.zeros(len(x))
for i in tid_data.index:
    start = tid_data.loc[i, 1]
    end = tid_data.loc[i, 2]
    for k in np.arange(start, end, 5):
        y_tid[k // 5] = 1
for i in epod_data.index:
    start = epod_data.loc[i, 3]
    end = epod_data.loc[i, 4]
    for k in np.arange(start, end, 5):
        y_epod[k // 5] = 1
epod_only = 0
tid_only = 0
both = 0
none = 0
for i in x:
    k = i // 5
    if (y_epod[k] == 1) & (y_tid[k] == 1):
        both += 5
    elif (y_epod[k] == 0) & (y_tid[k] == 1):
        tid_only += 5
    elif (y_epod[k] == 1) & (y_tid[k] == 0):
        epod_only += 5
    elif (y_epod[k] == 0) & (y_tid[k] == 0):
        none += 5

venn2(
    subsets=(epod_only // 1000, tid_only // 1000, both // 1000),
    set_labels=("EPODs", "TIDs"),
)
plt.savefig(venn_plot)
plt.close()
