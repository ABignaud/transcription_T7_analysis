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
import bacchus.plot as bcp

# Import snakemake.
wt_file = snakemake.input.mat_wt
pompA_file = snakemake.input.mat_pompA
prpsM_file = snakemake.input.mat_prpsM
cmap = snakemake.params.cmap
res = snakemake.params.res
width = snakemake.params.width

out_wt = str(snakemake.output.wt)
out_pompA = str(snakemake.output.pompA)
out_prpsM = str(snakemake.output.prpsM)

# Make sure output directory exists.
os.makedirs(str((snakemake.params.outdir)), exist_ok=True)

files = [wt_file, pompA_file, prpsM_file]
outs = [out_wt, out_pompA, out_prpsM]
title = ['WT', "pompA", "prpsM"]
hic = [0, 0, 0]
start = 300000 // res
end = 420000 // res
subsample = 24593739
for i, file in enumerate(files):
    clr = cooler.Cooler(f"{file}::/resolutions/{res}")
    mat = clr.matrix(balance=False, sparse=True)[: :]
    print(mat.sum())
    mat = hcs.subsample_contacts(mat, subsample) 
    mat = hcs.normalize_sparse(mat, norm='ICE', iterations=100, n_mad=10)
    hic[i] = mat.toarray()
    bcp.contact_map(
        hic[i],
        binning=res,
        title=title[i],
        zmax=0.002,
        out_file=outs[i],
        start=start,
        end=end,
    )
