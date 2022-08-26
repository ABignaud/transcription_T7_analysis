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

wildcard_constraints:
  hic_library = "|".join(samples.library[samples.type == 'hic']),
  rna_library = "|".join(samples.library[samples.type == 'rna']),
  rna_se_library = "|".join(samples.library[(samples.type == 'rna') & (samples.sequencing == 'SE')]),
  species = "|".join(species),
  condition = "|".join(conditions),

# Pipeline sub-workflows
include: 'rules/00_common.smk'
include: 'rules/01_hic_processing.smk'
include: 'rules/02_rna_processing.smk'
# include: 'rules/03_gapr_processing.smk'
include: 'rules/04_ecoli_matrix_analysis.smk'
# include: 'rules/04_pileup_analysis.smk'

rule all:
    input:
        # 01
        expand(join(OUT_DIR, 'HiC', '{condition}.mcool'), condition=conditions),
        # 02
        expand(join(OUT_DIR, 'RNA_tracks', '{rna_library}_unstranded.bw'), rna_library=samples[samples.type == "rna"].index),
        expand(join(OUT_DIR, 'RNA_tracks', '{condition}_rna_unstranded.bw'), condition=conditions),
        # 03
        # 04
        expand(join(OUT_DIR, 'figures', 'E_coli_contact_map', 'WT_{res}kb.pdf'), res=[1, 2, 4]),
        expand(
          join(OUT_DIR, 'figures', 'E_coli_contact_map', 'mat_zoom', 'region_{pos1}_{pos2}.pdf'),
          zip,
          **{'pos1': [420, 1110, 1960, 3436, 3870], 'pos2': [484, 1174, 2024, 3500, 3934]},
        ),
        expand(join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_map.pdf'), species=species),
        expand(join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_corr.pdf'), species=species),
        expand(join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_pileup.pdf'), species=species),
        expand(join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_pileup_pos_TU.pdf'), species=species),
        expand(join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_pileup_neg_TU.pdf'), species=species),
        join(OUT_DIR, 'figures', 'E_coli_contact_map', 'CIDs', 'comaprison_Lioy.pdf'),

        
        

