#!/bin/env snakemake -s

# Rules to generates contact map plot.

# Generate GC content using dnaglider
rule generate_gc:
    input:
        fasta = join(REF_DIR, config['E_coli']['ref'])
    output:
        gc = join(OUT_DIR, 'E_coli', 'ref', 'gc_content.tsv')
    threads: 1
    conda: "../envs/split_fastq.yaml"
    shell: 
        """
        seqkit sliding \
            --step 100 \
            --window 500 \
            --threads {threads} \
            {input.fasta} | \
                seqkit fx2tab --gc | \
                cut -f1,4 > {output.gc}
        """


# Rule plot the matrix.
rule plot_matrix:
    input:
        mat = join(OUT_DIR, 'HiC', 'E_coli.mcool'),
    params:
        res = lambda w: int(w.res) * 1000,
        cmap = config['cmap'],
        title = 'E. coli WT - binning {res}kb',
        outdir = join(OUT_DIR, 'figures', 'WT_contact_map'),
    output:
        join(OUT_DIR, 'figures', 'WT_contact_map', 'WT_{res}kb.pdf'),
    conda: '../envs/bacchus.yaml'
    threads: 1
    script: '../scripts/plot_matrix.py'


# EPODs positions download from Freddolino et al., PLOS, 2021.
# Rule plot the zoom of matrices with different metrics.
rule plot_zoom:
    input:
        annotation = join(REF_DIR, config['E_coli']['annotation']),
        cool_wt = join(OUT_DIR, 'HiC', 'E_coli.mcool'),
        cool_rf = join(OUT_DIR, 'HiC', 'E_coli_rif.mcool'),
        rna_wt = join(OUT_DIR, 'RNA_tracks', 'E_coli_rna_unstranded.bw'),
        rna_rf = join(OUT_DIR, 'RNA_tracks', 'E_coli_rif_rna_unstranded.bw'),
        cov_wt = join(OUT_DIR, 'HiC_tracks', 'CC328.bw'),
        cov_rf = join(OUT_DIR, 'HiC_tracks', 'CC419.bw'),
        gc = join(OUT_DIR, 'E_coli', 'ref', 'gc_content.tsv'),
        EPODs = join(REF_DIR, config['epods']),
    params:
        res = config['E_coli']['res'][0],
        cmap = config['cmap'],
        width = config['zoom_width'],
        positions = lambda w: [w.pos1, w.pos2],
        outdir = join(OUT_DIR, 'figures', 'WT_contact_map', 'mat_zoom'),
    output:
        join(OUT_DIR, 'figures', 'WT_contact_map', 'mat_zoom', 'region_{pos1}_{pos2}.pdf'),
    conda: '../envs/bacchus.yaml'
    threads: 1
    script: '../scripts/plot_zoom.py'


# Run Macrodomain and CIDs detection using directional index and comapre it to 
# Lioy dataset.
rule CIDs_comparison_Lioy:
    input:
        mat = join(OUT_DIR, 'HiC', 'E_coli.mcool'),
    params:
        res = 5000,
        cmap = config['cmap'],
        outdir = join(OUT_DIR, 'figures', 'WT_contact_map', 'CIDs'),
    output:
        macrodomain = join(OUT_DIR, 'figures', 'WT_contact_map', 'CIDs', 'macrodomaim.pdf'),
        cid = join(OUT_DIR, 'figures', 'WT_contact_map', 'CIDs', 'CIDs.pdf'),
        compare = join(OUT_DIR, 'figures', 'WT_contact_map', 'CIDs', 'comparison_Lioy.pdf'),
    conda: '../envs/bacchus.yaml'
    threads: 1
    script: '../scripts/CIDs_analysis_Lioy.py'


# Run high resolution CIDs analysis based on Directional index and insulation 
# score.
rule CIDs_analysis:
    input:
        mat = join(OUT_DIR, 'HiC', 'E_coli.mcool'),
        rna = join(OUT_DIR, 'RNA_tracks', 'E_coli_rna_unstranded.bw'),
        annotation = join(REF_DIR, config['E_coli']['annotation']),
    params:
        cmap = config['cmap'],
        dpi = 100,
        outdir = join(OUT_DIR, 'figures', 'WT_contact_map', 'CIDs'),
    output:
        full_mat = join(OUT_DIR, 'figures', 'WT_contact_map', 'CIDs', 'full_mat_CIDs.pdf'),
        zoom1 = join(OUT_DIR, 'figures', 'WT_contact_map', 'CIDs', 'zoom1_CIDs.pdf'),
        zoom2 = join(OUT_DIR, 'figures', 'WT_contact_map', 'CIDs', 'zoom2_CIDs.pdf'),
        full_mat_di = join(OUT_DIR, 'figures', 'WT_contact_map', 'CIDs', 'full_mat_CIDs_DI.pdf'),
        zoom1_di = join(OUT_DIR, 'figures', 'WT_contact_map', 'CIDs', 'zoom1_CIDs_DI.pdf'),
        zoom2_di = join(OUT_DIR, 'figures', 'WT_contact_map', 'CIDs', 'zoom2_CIDs_DI.pdf'),
        bor_trans = join(OUT_DIR, 'figures', 'WT_contact_map', 'CIDs', 'genes_transcription_borders.pdf'),
        trans_bor = join(OUT_DIR, 'figures', 'WT_contact_map', 'CIDs', 'genes_transcription_borders_distance_mean.pdf'),
        text_file = join(OUT_DIR, 'figures', 'WT_contact_map', 'CIDs', 'CIDs_statistics.txt'),
    conda: '../envs/bacchus.yaml'
    threads: 1
    script: '../scripts/CIDs_analysis_resolution.py'


# Detect TIDs and compute some metrics on them 
rule TIDs_analysis:
    input:
        mat = join(OUT_DIR, 'HiC', 'E_coli.mcool'),
        rna = join(OUT_DIR, 'RNA_tracks', 'E_coli_rna_unstranded.bw'),
        cov = join(OUT_DIR, 'HiC_tracks', 'CC328.bw'),
        gc = join(OUT_DIR, 'E_coli', 'ref', 'gc_content.tsv'),
        frags = join(TMP, 'hicstuff', 'CC328', 'CC328.frags.tsv'),
        EPODs = join(REF_DIR, config['epods']),
    params: outdir = join(OUT_DIR, 'figures', 'WT_contact_map', 'TIDs'),
    output:
        bed = join(OUT_DIR, 'figures', 'WT_contact_map', 'TIDs', 'TIDs.bed'),
        text_file = join(OUT_DIR, 'figures', 'WT_contact_map', 'TIDs', 'TIDs_statistics.txt'),
        gc = join(OUT_DIR, 'figures', 'WT_contact_map', 'TIDs', 'TIDs_GC.pdf'),
        cov = join(OUT_DIR, 'figures', 'WT_contact_map', 'TIDs', 'TIDs_cov.pdf'),
        RS = join(OUT_DIR, 'figures', 'WT_contact_map', 'TIDs', 'TIDs_RS.pdf'),
        rna = join(OUT_DIR, 'figures', 'WT_contact_map', 'TIDs', 'TIDs_RNA.pdf'),
        dist_bar = join(OUT_DIR, 'figures', 'WT_contact_map', 'TIDs', 'TIDs_bar_dist.pdf'),
        dist_line = join(OUT_DIR, 'figures', 'WT_contact_map', 'TIDs', 'TIDs_line_dist.pdf'),
        venn_plot = join(OUT_DIR, 'figures', 'WT_contact_map', 'TIDs', 'TIDs_Venn.pdf'),
    conda: '../envs/bacchus.yaml'
    threads: 1
    script: '../scripts/TIDs_analysis.py'

rule all_fig1:
    input:
        expand(join(OUT_DIR, 'figures', 'WT_contact_map', 'WT_{res}kb.pdf'), res=[1, 2, 5]),
        expand(
          join(OUT_DIR, 'figures', 'WT_contact_map', 'mat_zoom', 'region_{pos1}_{pos2}.pdf'),
          zip,
          **{'pos1': [420, 1110, 1960, 3436, 3870], 'pos2': [484, 1174, 2024, 3500, 3934]},
        ),
        # join(OUT_DIR, 'figures', 'WT_contact_map', 'CIDs', 'comparison_Lioy.pdf'),
        # join(OUT_DIR, 'figures', 'WT_contact_map', 'CIDs', 'full_mat_CIDs.pdf'),
        join(OUT_DIR, 'figures', 'WT_contact_map', 'TIDs', 'TIDs_statistics.txt'),
    output: touch(join(TMP, 'all_fig1.done'))
    threads: 1
