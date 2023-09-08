#!/bin/env snakemake -s

# Rules to plot matrices and ratio of the contact map of mutants or drug treated
# bacteria.

# Rule plot the matrix topA.
rule plot_matrix_topA:
    input:
        mat_ara = join(OUT_DIR, 'HiC', 'T7_ara.mcool'),
        mat_topA = join(OUT_DIR, 'HiC', 'T7_topA.mcool'),
        mat_ara_rif = join(OUT_DIR, 'HiC', 'T7_ara_rif.mcool'),
        mat_topA_rif = join(OUT_DIR, 'HiC', 'T7_topA_rif.mcool'),
    params:
        cmap = config['cmap'],
        res = 1000,
        res_ratio = 2000,
        width = 50000,
        outdir = join(OUT_DIR, 'figures', 'T7_system', 'TopA'),
    output:
        matrices = join(OUT_DIR, 'figures', 'T7_system', 'TopA', 'matrices.pdf'),
        ratio = join(OUT_DIR, 'figures', 'T7_system', 'TopA', 'ratio.pdf'),
        ratio_rif = join(OUT_DIR, 'figures', 'T7_system', 'TopA', 'ratio_rif.pdf'),
        # zoom = join(OUT_DIR, 'figures', 'T7_system', 'TopA', 'zoom_matrices.pdf'),
        zoom_ratio = join(OUT_DIR, 'figures', 'T7_system', 'TopA', 'zoom_ratio.pdf'),
    conda: '../envs/bacchus.yaml'
    threads: config['threads_serp']
    script: '../scripts/topA_matrices.py'

# Rule plot the matrix novobiocin.
rule plot_matrix_novo:
    input:
        mat_ara = join(OUT_DIR, 'HiC', 'T7_ara.mcool'),
        mat_novo = join(OUT_DIR, 'HiC', 'T7_novo.mcool'),
        mat_ara_rif = join(OUT_DIR, 'HiC', 'T7_ara_rif.mcool'),
        mat_novo_rif = join(OUT_DIR, 'HiC', 'T7_novo_rif.mcool'),
    params:
        cmap = config['cmap'],
        res = 1000,
        res_ratio = 2000,
        width = 50000,
        outdir = join(OUT_DIR, 'figures', 'T7_system', 'novo'),
    output:
        matrices = join(OUT_DIR, 'figures', 'T7_system', 'novo', 'matrices.pdf'),
        ratio = join(OUT_DIR, 'figures', 'T7_system', 'novo', 'ratio.pdf'),
        ratio_rif = join(OUT_DIR, 'figures', 'T7_system', 'novo', 'ratio_rif.pdf'),
        # zoom = join(OUT_DIR, 'figures', 'T7_system', 'novo', 'zoom_matrices.pdf'),
        zoom_ratio = join(OUT_DIR, 'figures', 'T7_system', 'novo', 'zoom_ratio.pdf'),
    conda: '../envs/bacchus.yaml'
    threads: config['threads_serp']
    script: '../scripts/novobiocin_matrices.py'

# Do the correlations between hic and rna tracks for ara induction.
rule corr_T7_ara:
    input:
        hic = join(OUT_DIR, 'HiC', 'T7_ara.mcool'),
        rna = join(OUT_DIR,  'RNA_tracks', 'T7_ara_rna_unstranded.bw'),
    output:
        mat_fig = join(OUT_DIR, 'figures', 'T7_system', 'TopA', 'ara_map.pdf'),
        corr_fig = join(OUT_DIR, 'figures', 'T7_system', 'TopA', 'rna_ara_corr.pdf'),
    params:
        binning = 500,
        species = 'T7_ara',
        circular = True,
    conda: '../envs/bacchus.yaml'
    threads: 1
    script: '../scripts/corr_hic_rna.py'

# Generates pileups plot of HEGs from the output of previous tracks for each 
# topA mutant and ara condition.
rule pileup_topA:
    input:
        hic = join(OUT_DIR, 'HiC', '{condition}.mcool'),
        rna = join(OUT_DIR,  'RNA_tracks', 'T7_ara_rna_unstranded.bw'),
        annotation = join(REF_DIR, 'RSG_B834', 'RSG_B834.gff'),
    output:
        pileup_pos_TU = join(OUT_DIR, 'figures', 'T7_system', 'TopA', '{condition}', '{condition}_{res}_pileup_pos_TU.pdf'),
        pileup_neg_TU = join(OUT_DIR, 'figures', 'T7_system', 'TopA', '{condition}', '{condition}_{res}_pileup_neg_TU.pdf'),
        pileup = join(OUT_DIR, 'figures', 'T7_system', 'TopA', '{condition}', '{condition}_{res}_pileup.pdf'),
        pileup_zoom = join(OUT_DIR, 'figures', 'T7_system', 'TopA', '{condition}', '{condition}_{res}_pileup_zoom.pdf'),
        pileup_zoom2 = join(OUT_DIR, 'figures', 'T7_system', 'TopA', '{condition}', '{condition}_{res}_pileup_zoom2.pdf'),
        rpkm = join(OUT_DIR, 'figures', 'T7_system', 'TopA', '{condition}', '{condition}_{res}_rpkm.pdf'),
    params:
        binning = '{res}',
        species = '{condition}',
        out_dir = join(OUT_DIR, 'figures', 'T7_system', 'TopA', '{condition}'),
        circular = lambda w: config["E_coli"]['circular'],
        window = lambda w: config["E_coli"]['window'],
        threshold = lambda w: config["E_coli"]['threshold'],
        unit_length = lambda w: config["E_coli"]['unit_length'],
    conda: '../envs/bacchus.yaml'
    threads: 1
    script: '../scripts/pileup_species.py'

# Endogenous promoters:
rule endogenous_promoters:
    input:
        mat_wt = join(OUT_DIR, 'HiC', 'E_coli.mcool'),
        mat_pompA = join(OUT_DIR, 'HiC', 'pompA.mcool'),
        mat_prpsM = join(OUT_DIR, 'HiC', 'prpsM.mcool'),
    params:
        cmap = config['cmap'],
        res = 1000,
        width = 32000,
        outdir = join(OUT_DIR, 'figures', 'T7_system', 'prom'),
    output:
        wt = join(OUT_DIR, 'figures', 'T7_system', 'prom', 'wt.pdf'), 
        pompA = join(OUT_DIR, 'figures', 'T7_system', 'prom', 'pompA.pdf'),
        prspM = join(OUT_DIR, 'figures', 'T7_system', 'prom', 'prpsM.pdf'),
    conda: '../envs/bacchus.yaml'
    threads: 1
    script: '../scripts/prom_matrices.py'


rule all_drugs:
    input:
        matrices = join(OUT_DIR, 'figures', 'T7_system', 'TopA', 'matrices.pdf'),
        zoom = join(OUT_DIR, 'figures', 'T7_system', 'TopA', 'zoom_matrices.pdf'),
        zoom_ratio = join(OUT_DIR, 'figures', 'T7_system', 'TopA', 'zoom_ratio.pdf'),
        mat_fig = join(OUT_DIR, 'figures', 'T7_system', 'TopA', 'ara_map.pdf'),
        corr_fig = join(OUT_DIR, 'figures', 'T7_system', 'TopA', 'rna_ara_corr.pdf'),
        matrices_novo = join(OUT_DIR, 'figures', 'T7_system', 'novo', 'matrices.pdf'),
        ratio_novo = join(OUT_DIR, 'figures', 'T7_system', 'novo', 'ratio.pdf'),
        ratio_rif_novo = join(OUT_DIR, 'figures', 'T7_system', 'novo', 'ratio_rif.pdf'),
        zoom_ratio_novo = join(OUT_DIR, 'figures', 'T7_system', 'novo', 'zoom_ratio.pdf'),
        pileup_ara = join(OUT_DIR, 'figures', 'T7_system', 'TopA', 'T7_ara', 'T7_ara_1000_pileup.pdf'),
        pileup_topA = join(OUT_DIR, 'figures', 'T7_system', 'TopA', 'T7_topA', 'T7_topA_1000_pileup.pdf'),
        pileup_ara5 = join(OUT_DIR, 'figures', 'T7_system', 'TopA', 'T7_ara', 'T7_ara_500_pileup.pdf'),
        pileup_topA5 = join(OUT_DIR, 'figures', 'T7_system', 'TopA', 'T7_topA', 'T7_topA_500_pileup.pdf'),wt = join(OUT_DIR, 'figures', 'T7_system', 'prom', 'wt.pdf'), 
        pompA = join(OUT_DIR, 'figures', 'T7_system', 'prom', 'pompA.pdf'),
        prspM = join(OUT_DIR, 'figures', 'T7_system', 'prom', 'prpsM.pdf'),
    output:
        touch(join(TMP, 'drugs_fig.done'))
    threads: 1


# rule chlorophenicol:
# rule stop_codon:
# rule matP:

