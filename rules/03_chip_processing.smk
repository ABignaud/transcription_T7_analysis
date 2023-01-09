#!/bin/env snakemake -s

# Rule to generates the ChIP tracks. Based on timymapper pipeline. 
# https://github.com/js2264/tinyMapper.git

rule align_chip:
  input:
    index_flag = lambda w: join(TMP, 'ref', f'{samples.species[w.chip_library]}_bt2_index.done'),
    R1 = join(TMP, 'split_reads', '{chip_library}_R1', '{chip_library}_R1.{split}.fq.gz'),
    R2 = join(TMP, 'split_reads', '{chip_library}_R2', '{chip_library}_R2.{split}.fq.gz'),
  output: 
    join(TMP, 'align', '{chip_library}', '{chip_library}.{split}.bam'),
  params:
    index = lambda w: join(TMP, 'ref', f'{samples.species[w.chip_library]}_genome'),
    bt2_presets = config['bowtie2_args'],
  threads: config['threads_small']
  conda: "../envs/gen_tracks.yaml"
  log: join(OUT_DIR, 'logs', 'rna', '{chip_library}_{split}.log')
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

rule merge_split_chip_alignments:
  input:
    expand(
      join(TMP, 'align', '{{chip_library}}',  '{{chip_library}}.{split}.bam'),
      split=split_names
    )
  output: join(TMP, 'bam', '{chip_library}.bam')
  threads: config['threads_small']
  conda: "../envs/gen_tracks.yaml"
  shell:
    """
    samtools merge -n -O BAM -@ {threads} - {input} | 
        samtools sort -@ {threads} --output-fmt bam -l 9 -o {output}
    """

rule chip_coverage:
    input:
        bam = join(TMP, 'bam', '{chip_library}.bam'),
    output:
        unstranded = join(OUT_DIR, 'ChIP_tracks', '{chip_library}.bw'),
    threads: config['threads_small']
    conda: "../envs/gen_tracks.yaml"
    shell:
        """
        samtools index {input.bam} -@ {threads}
        bamCoverage --bam {input.bam} \
            --outFileName {output.unstranded} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads 
        """
        