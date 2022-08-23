#!/bin/env snakemake -s

# Rule to generates the RNA tracks


# Map the reads and sort them.
rule align_rna:
  input:
    index_flag = join(OUT_DIR, '{species}', 'ref', 'bt2_index.done'),
    R1 = join(TMP, 'split_reads', '{rna_library}_R1', '{rna_library}_R1.{split}.fq.gz'),
    R2 = join(TMP, 'split_reads', '{rna_library}_R2', '{rna_library}_R2.{split}.fq.gz'),
  output: 
    join(TMP, '{species}', 'split_reads', 'rna_bam',  '{rna_library}.{split}.bam'),
  params:
    index = lambda w: join(OUT_DIR, samples.species[w.rna_library], 'ref' , 'genome'),
    bt2_presets = config['bowtie2_args'],
  threads: config['threads']
#   singularity: "docker://koszullab/hicstuff:v3.1.0"
  conda: "../envs/gen_tracks.yaml"
  log: join(OUT_DIR, '{species}', 'logs', 'mapping', '{rna_library}_{split}.log')
  shell:
    """
    bowtie2 {params.bt2_presets} \
            -p {threads} \
            -x {params.index} \
            --maxins 1000 \
            -1 {input.R1} \
            -2 {input.R2} 2> {log} | \
        samtools sort -@ {threads} -n -O BAM - | \
        samtools fixmate -@ {threads} --output-fmt bam -m - - | \
        samtools sort -@ {threads} -O BAM - | \
        samtools markdup -@ {threads} --output-fmt bam -r - - | \
        samtools view -@ {threads} --output-fmt bam -f 2 -q 10 -1 -b - | \
        samtools sort -@ {threads} --output-fmt bam -l 9 -o {output}
    """

rule align_SE_rna:
  input:
    index_flag = join(OUT_DIR, '{species}', 'ref', 'bt2_index.done'),
    fq = join(TMP, 'split_reads', '{rna_se_library}', '{rna_se_library}.{split}.fq.gz'),
  output: 
    join(TMP, '{species}', 'split_reads', 'rna_bam',  '{rna_se_library}.{split}.bam'),
  params:
    index = lambda w: join(OUT_DIR, samples.species[w.rna_se_library], 'ref' , 'genome'),
    bt2_presets = config['bowtie2_args'],
  threads: config['threads']
#   singularity: "docker://koszullab/hicstuff:v3.1.0"
  conda: "../envs/gen_tracks.yaml"
  log: join(OUT_DIR, '{species}', 'logs', 'mapping', '{rna_se_library}_{split}.log')
  shell:
    """
    bowtie2 {params.bt2_presets} \
            -p {threads} \
            -x {params.index} \
            --maxins 1000 \
            -U {input.fq} 2> {log} | \
        samtools sort -@ {threads} -O BAM - | \
        samtools markdup -@ {threads} --output-fmt bam -r - - | \
        samtools view -@ {threads} --output-fmt bam -q 10 -1 -b - | \
        samtools sort -@ {threads} --output-fmt bam -l 9 -o {output}
    """

rule merge_split_rna_alignments:
  input:
    expand(
      join(TMP, '{{species}}',  'split_reads', 'rna_bam', '{{rna_library}}.{split}.bam'),
      split=split_names
    )
  output: join(OUT_DIR, '{species}', 'bam', '{rna_library}.bam')
  threads: config['threads']
  # singularity: "docker://biocontainers/samtools:v1.7.0_cv4"
  conda: "../envs/gen_tracks.yaml"
  shell:
    """
    samtools merge -n -O BAM -@ {threads} - {input} | 
        samtools sort -@ {threads} --output-fmt bam -l 9 -o {output}
    """


# Build RNA tracks
rule rna_coverage:
    input:
        bam = join(OUT_DIR, '{species}', 'bam', '{rna_library}.bam'),
    output:
        unstranded = join(OUT_DIR, '{species}',  'tracks', '{rna_library}_unstranded.bw'),
        fw = join(OUT_DIR, '{species}', 'tracks', '{rna_library}_forward.bw'),
        rv = join(OUT_DIR, '{species}', 'tracks', '{rna_library}_reverse.bw'),
    threads: config['threads']
    conda: "../envs/gen_tracks.yaml"
    shell:
        """
        samtools index {input.bam} -@ {threads}
        bamCoverage --bam {input.bam} \
            --outFileName {output.unstranded} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads 200 \
            --ignoreDuplicates
        bamCoverage --bam {input.bam} \
            --outFileName {output.fw} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads 200 \
            --ignoreDuplicates\
            --filterRNAstrand forward
        bamCoverage --bam {input.bam} \
            --outFileName {output.rv} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads 200 \
            --ignoreDuplicates\
            --filterRNAstrand reverse
        """


# Merge multiple bam files into one.
rule merge_rna_bam:
    input:
        lambda w: [join(
                OUT_DIR, '{species}', 'bam', f'{i}.bam'
            ) for i in samples[
                (samples.species == w.species) & (samples.type == 'rna')
            ].index]
    output: join(OUT_DIR, '{species}', 'bam', 'rna_merge.bam')
    threads: config["threads"]
    conda: "../envs/gen_tracks.yaml"
    params: 
        lambda w: " ".join([join(
                OUT_DIR, w.species, 'bam', f'{i}.bam'
            ) for i in samples[
                (samples.species == w.species) & (samples.type == 'rna')
            ].index])
    shell:'samtools merge -l 9 -O BAM -@ {threads} -f {output} {params}'
        

# RNA coverage for species merge.
rule rna_coverage_merge:
    input:
        join(OUT_DIR, '{species}', 'bam', 'rna_merge.bam'),
    output:
        unstranded = join(OUT_DIR,  '{species}', 'tracks', '{species}_rna_unstranded.bw'),
        fw = join(OUT_DIR,  '{species}', 'tracks', '{species}_rna_forward.bw'),
        rv = join(OUT_DIR,  '{species}', 'tracks', '{species}_rna_reverse.bw'),
    threads: config['threads']
    conda: "../envs/gen_tracks.yaml"
    shell:
        """
        samtools index {input} -@ {threads}
        bamCoverage --bam {input} \
            --outFileName {output.unstranded} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads 200 \
            --ignoreDuplicates
        bamCoverage --bam {input} \
            --outFileName {output.fw} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads 200 \
            --ignoreDuplicates\
            --filterRNAstrand forward
        bamCoverage --bam {input} \
            --outFileName {output.rv} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads 200 \
            --ignoreDuplicates \
            --filterRNAstrand reverse
        """