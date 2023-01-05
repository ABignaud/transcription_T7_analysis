#!/bin/env snakemake -s

# snakemake --rulegraph | dot -Tsvg > images/rulegraph.svg
# snakemake --dag | dot -Tsvg > images/dag.svg

# This file can be run using snakemake. It was tested on snakemake 5.10.0.
# It orchestrates the analysis of the gene trasncription impact on the HiC map
# from different species.

import numpy as np
import pandas as pd
from os.path import join
from snakemake.utils import validate

# Set parameters.
shell.prefix("set -euo pipefail;")

# LOAD CONFIG FILES
configfile: 'config/config.yaml'

samples = pd.read_csv(
    config['samples'], 
    sep=';', 
    dtype=str,
    comment='#',
).set_index(['library'], drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

species = np.unique(samples.species)
conditions = np.unique(samples.condition)

OUT_DIR = join(config['base_dir'], config['out_dir'])
TMP = join(config['base_dir'], config['tmp_dir'])
REF_DIR = join(config['base_dir'], config['ref_dir'])
FASTQ_DIR = join(config['base_dir'], config['fastq_dir'])

CONTROL_CHIP = 'CCchip14'

wildcard_constraints:
    hic_library = "|".join(samples.library[samples.type == 'hic']),
    rna_library = "|".join(samples.library[samples.type == 'rna']),
    chip_library = f'{"|".join(samples.library[samples.type == "chip"])}|{CONTROL_CHIP}',
    rna_se_library = "|".join(samples.library[(samples.type == 'rna') & (samples.sequencing == 'SE')]),
    species = "|".join(species),
    condition = "|".join(conditions),

# Pipeline sub-workflows
include: 'rules/00_common.smk'
include: 'rules/01_hic_processing.smk'
include: 'rules/02_rna_processing.smk'
include: 'rules/03_chip_processing.smk'
include: 'rules/04_WT_ecoli_analysis.smk'
# include: 'rules/04_pileup_analysis.smk'

rule all:
    input:
        # 01 - HiC
        expand(join(OUT_DIR, 'HiC', '{condition}.mcool'), condition=conditions),
        # 02 - RNA
        expand(
            join(OUT_DIR, 'RNA_tracks', '{condition}_rna_unstranded.bw'),
            condition=conditions
        ),
        # 03 - ChIP
        expand(join(
                OUT_DIR, 'ChIP_tracks', '{chip_library}.bw'
            ),
            chip_library = samples.library[samples.type == 'chip'],
        ),
        # 04 - Fig1
        # join(TMP, 'all_fig1.done'),
        # 05 - Fig2
        # 06 - Fig3
        # 07 - Fig4 - pileup
        # join(TMP, 'pileup.done'),
        # 08 - Others

        
        

