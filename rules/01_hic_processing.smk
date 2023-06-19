#!/bin/env snakemake -s

# Rules to generates the HiC matrices.


# Alignment of a single fastq split from a Hi-C library with hicstuff iteralign.
rule split_align_hic:
  input:
      index_flag = lambda w: join(TMP, 'ref', f'{samples.species[w.hic_library]}_bt2_index.done'),
      fq = join(TMP, 'split_reads', '{hic_library}_R{end}', '{hic_library}_R{end}.{split}.fq.gz'),
  output: join(TMP, 'align', '{hic_library}_R{end}', '{hic_library}_R{end}.{split}.bam')
  params:
      index = lambda w: join(TMP, 'ref', f'{samples.species[w.hic_library]}_genome'),
      read_len = lambda w: samples.read_len[w.hic_library],
      tmp = join(TMP, 'split_reads', '{hic_library}_R{end}_{split}'),
  threads: config['threads_small']
  conda: "../envs/hic_processing.yaml"
  shell:
    """
    hicstuff iteralign \
        --aligner bowtie2 \
        --read-len {params.read_len} \
        --threads {threads} \
        --tempdir {params.tmp} \
        --genome {params.index} \
        --out-bam {output} \
        {input.fq}
    """


# Merge splits from individual mapping jobs into one bam file per library end
rule merge_split_hic_alignments:
  input:
      expand(
          join(TMP, 'align', '{{hic_library}}_R{{end}}', '{{hic_library}}_R{{end}}.{split}.bam'),
          split=split_names
      )
  output: join(TMP, 'bam', '{hic_library}_R{end}.bam')
  threads: config['threads']
  conda: "../envs/hic_processing.yaml"
  shell: "samtools merge -n -O BAM -@ {threads} {output} {input}"


rule hic_coverage:
    input:
        bam = join(TMP, 'bam', '{hic_library}_R1.bam'),
    params:
        sorted_bam = join(TMP, 'bam', '{hic_library}_R1_sorted.bam'),
    output:
        unstranded = join(OUT_DIR, 'HiC_tracks', '{hic_library}.bw'),
    threads: config['threads_small']
    conda: "../envs/gen_tracks.yaml"
    shell:
        """
        samtools sort {input.bam} -@ {threads} -o {params.sorted_bam}
        samtools index {params.sorted_bam} -@ {threads}
        bamCoverage --bam {params.sorted_bam} \
            --outFileName {output.unstranded} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM 
        """


## 00 Generate Hi-C pairs files
rule generate_pairs:
  input:
      bam1 = join(TMP, 'bam', "{hic_library}_R1.bam"),
      bam2 = join(TMP, 'bam', "{hic_library}_R2.bam"),
  output:
      cool = join(TMP, 'hicstuff', '{hic_library}', '{hic_library}.cool'),
      frags = join(TMP, 'hicstuff', '{hic_library}', '{hic_library}.frags.tsv'),
  params:
      enz = lambda w: samples.enzyme[w.hic_library],
      hicdir = directory(join(TMP, 'hicstuff', '{hic_library}')),
      idx = lambda w: join(TMP, 'ref', f'{samples.species[w.hic_library]}_genome'),
      log = lambda w: join(OUT_DIR, samples.species[w.hic_library], f"logs/hicstuff/{w.hic_library}.log"),
      log_dir = lambda w: join(OUT_DIR, samples.species[w.hic_library], "logs/hicstuff/"),
  threads: 1
  conda: "../envs/hic_processing.yaml"
  shell:
    """
    mkdir -p {params.log_dir}
    hicstuff pipeline --force \
                      -e {params.enz} \
                      -g {params.idx} \
                      -o {params.hicdir} \
                      -S bam \
                      -P {wildcards.hic_library} \
                      -npD \
                      -M cool \
                      {input.bam1} \
                      {input.bam2} 2> {params.log}
    """

rule merge_rebin_cool:
  input:
      lambda w: [join(TMP, 'hicstuff', f'{i}', f'{i}.cool') for i in samples[(samples.condition == w.condition) & (samples.type == 'hic')].index]
  params:
      input = lambda w: " ".join([join(TMP, 'hicstuff', f'{i}', f'{i}.cool') for i in samples[(samples.condition == w.condition) & (samples.type == 'hic')].index]),
      prefix = join(TMP, 'cool', '{condition}_500bp'),
  output: join(TMP, 'cool', '{condition}_500bp.cool')
  threads: 1
  conda: "../envs/hic_processing.yaml"
  shell: 
      """
      hicstuff rebin -b 500bp {params.input} {params.prefix} 
      cooler merge {params.output} {params.prefix}.cool
      """


# Generate multiple resolutions and normalize mcool files.
rule zoomify_normalize_cool:
  input:
      cool = join(TMP, 'cool', '{condition}_500bp.cool')
  output: join(OUT_DIR, 'HiC', '{condition}.mcool')
  threads: config['threads_small']
  params:
      res = lambda w: ",".join(map(str, config[samples.species[samples.condition == w.condition][0]]['res'])),
      balance_args = config['balance_args'],
  conda: "../envs/hic_processing.yaml"
  shell:
      """
      cooler zoomify -r {params.res} \
          --balance \
          --balance-args '{params.balance_args}' \
          -n {threads} \
          -o {output} \
          {input.cool}
      """
