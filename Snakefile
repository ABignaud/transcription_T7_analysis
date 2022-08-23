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

OUT_DIR = join(config['base_dir'], config['out_dir'])
TMP = join(config['base_dir'], config['tmp_dir'])
REF_DIR = join(config['base_dir'], config['ref_dir'])
FASTQ_DIR = join(config['base_dir'], config['fastq_dir'])

wildcard_constraints:
  hic_library = "|".join(samples.library[samples.type == 'hic']),
  rna_library = "|".join(samples.library[samples.type == 'rna']),
  rna_se_library = "|".join(samples.library[(samples.type == 'rna') & (samples.sequencing == 'SE')]),
  species = "|".join(species)

# Pipeline sub-workflows
include: 'rules/00_common.smk'
include: 'rules/01_hic_processing.smk'
include: 'rules/02_rna_processing.smk'
# include: 'rules/03_gapr_processing.smk'
include: 'rules/04_pileup_analysis.smk'

rule all:
    input:
        # 01
        expand(join(OUT_DIR, '{species}', 'cool', '{species}.mcool'), species=species),
        # 02
        expand(join(OUT_DIR, '{species}', 'tracks', '{rna_library}_unstranded.bw'), zip, **samples.species[samples.type == "rna"]),
        expand(join(OUT_DIR, '{species}', 'tracks', '{rna_library}_forward.bw'), zip, **samples.species[samples.type == "rna"]),
        expand(join(OUT_DIR, '{species}', 'tracks', '{rna_library}_reverse.bw'), zip, **samples.species[samples.type == "rna"]),
        expand(join(OUT_DIR, '{species}', 'tracks', '{species}_rna_unstranded.bw'), species=species),
        expand(join(OUT_DIR, '{species}', 'tracks', '{species}_rna_forward.bw'), species=species),
        expand(join(OUT_DIR, '{species}', 'tracks', '{species}_rna_reverse.bw'), species=species),
        # 03
        # 04
        # expand(join(OUT_DIR, '{species}', 'figures', '{species}_map.pdf'), species=species),
        # expand(join(OUT_DIR, '{species}', 'figures', '{species}_corr.pdf'), species=species),
        expand(join(OUT_DIR, '{species}', 'figures', '{species}_pileup.pdf'), species=species),
        expand(join(OUT_DIR, '{species}', 'figures', '{species}_pileup_pos_TU.pdf'), species=species),
        expand(join(OUT_DIR, '{species}', 'figures', '{species}_pileup_neg_TU.pdf'), species=species)

        
        

