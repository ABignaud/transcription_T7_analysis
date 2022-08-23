#!/bin/env snakemake -s

# Rules to generates the HiC matrices.


# Alignment of a single fastq split from a Hi-C library with hicstuff iteralign.
rule split_align_hic:
  input:
    index_flag = lambda w: join(OUT_DIR, samples.species[w.hic_library], 'ref', 'bt2_index.done'),
    fq = join(TMP, 'split_reads', '{hic_library}_R{end}', '{hic_library}_R{end}.{split}.fq.gz'),
  output: join(TMP, 'split_reads', '{hic_library}_R{end}', '{hic_library}_R{end}.{split}.bam')
  params:
    index = lambda w: join(OUT_DIR, samples.species[w.hic_library], 'ref', 'genome'),
    read_len = lambda w: samples.read_len[w.hic_library],
    tmp = join(TMP, 'split_reads', '{hic_library}_R{end}_{split}'),
    # log = lambda w: join(OUT_DIR, samples.species[w.hic_library], "logs/mapping/{hic_library}_{split}_{end}.log"),
  threads: config['threads']
  # singularity: "docker://koszullab/hicstuff:v3.1.0"
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
      join(TMP, 'split_reads', '{{hic_library}}_R{{end}}', '{{hic_library}}_R{{end}}.{split}.bam'),
      split=split_names
    )
  output: join(TMP, 'bam', '{hic_library}_R{end}.bam')
  threads: config['threads']
  # singularity: "docker://biocontainers/samtools:v1.7.0_cv4"
  conda: "../envs/hic_processing.yaml"
  shell: "samtools merge -n -O BAM -@ {threads} {output} {input}"


## 00 Generate Hi-C pairs files
rule generate_pairs:
  input:
    bam1 = join(TMP, 'bam', "{hic_library}_R1.bam"),
    bam2 = join(TMP, 'bam', "{hic_library}_R2.bam"),
  output:
    hicdir = directory(join(TMP, 'hicstuff', '{species}_{hic_library}')),
    pairs = join(OUT_DIR, '{species}', 'pairs', '{hic_library}.pairs'),
  params:
    enz = lambda w: samples.enzyme[w.hic_library],
    idx = lambda w: join(OUT_DIR, samples.species[w.hic_library], 'ref', 'genome'),
    log = lambda w: join(OUT_DIR, samples.species[w.hic_library], f"logs/hicstuff/{w.hic_library}.log"),
    log_dir = lambda w: join(OUT_DIR, samples.species[w.hic_library], "logs/hicstuff/"),
  threads: 1
  # singularity: "docker://koszullab/hicstuff:v3.1.0"
  conda: "../envs/hic_processing.yaml"
  shell:
    """
    mkdir -p {params.log_dir}
    hicstuff pipeline --force \
                      -e {params.enz} \
                      -g {params.idx} \
                      -o {output.hicdir} \
                      -S bam \
                      --filter \
                      -P {wildcards.hic_library} \
                      -npD \
                      {input.bam1} \
                      {input.bam2} 2> {params.log}
    cp {output.hicdir}/tmp/{wildcards.hic_library}.valid_idx_pcrfree.pairs {output.pairs}
    """


# Generates chrom_size file
rule chromsizes:
  input:  lambda w: join(REF_DIR, config[w.species]['ref'])
  output: join(OUT_DIR, '{species}', 'ref', 'chrom.sizes')
  conda: "../envs/hic_processing.yaml"
  threads: 1
  shell:
    """
    samtools faidx {input}
    cut -f1,2 "{input}.fai" > {output}
    """


# Convert pairs to cool
rule pairs_to_cool:
  input:
    chroms = lambda w: join(OUT_DIR, samples.species[w.hic_library], 'ref', 'chrom.sizes'),
    pairs = join(OUT_DIR, '{species}', 'pairs', '{hic_library}.pairs')
  output: join(OUT_DIR, '{species}', 'cool', '{hic_library}.cool')
  params:
    res = lambda w: config[w.species]['res'][0]
  conda: "../envs/hic_processing.yaml"
  shell:
    """
    cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \
      {input.chroms}:{params.res} \
      {input.pairs} \
      {output}                  
    """


rule merge_cool:
  input:
    lambda w: [join(OUT_DIR, w.species, 'cool', f'{i}.cool') for i in samples[(samples.species == w.species) & (samples.type == 'hic')].index]
  params:
    input = lambda w: " ".join([join(OUT_DIR, w.species, 'cool', f'{i}.cool') for i in samples[(samples.species == w.species) & (samples.type == 'hic')].index])
  output: join(TMP, 'cool', '{species}.cool')
  conda: "../envs/hic_processing.yaml"
  shell: 'cooler merge {output} {params.input}'


# Generate multiple resolutions and normalize mcool files.
rule zoomify_normalize_cool:
  input:
    cool = join(TMP, 'cool', '{species}.cool')
  output: join(OUT_DIR, '{species}', 'cool', '{species}.mcool')
  threads: config['threads']
  params:
    res = lambda w: ",".join(map(str, config[w.species]['res'])),
    balance_args = config['balance_args'],
  # singularity: "docker://cmdoret/cooler:0.8.5"
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

# Merge ends of Hi-C bam files and sort by coord
# rule merge_sort_bam_ends:
#   input: expand(join(TMP, 'bam', '{{hic_library}}_hic.{end}.bam'), end=['end1', 'end2'])
#   output: temporary(join(TMP, 'bam', '{hic_library}_hic.merged.bam'))
#   threads: NCPUS
#   conda: "../envs/hic_processing.yaml"
#   shell:
#     """
# 	samtools merge -n -@ {threads} -O BAM - {input} \
# 	  | samtools sort -@ {threads} -O BAM -o {output}
#     samtools index -@ {threads} {output}
#     """


# Visualise Hi-C coverage along genome for each library
# rule plot_hic_coverage:
#   input: join(TMP, 'bam', '{hic_library}_hic.end1.bam')
#   output:
#       plot = join(OUT, 'plots', 'coverage_hic_{hic_library}.pdf'),
# 	  text = join(OUT, 'cov_hic', 'coverage_hic_{hic_library}.bedgraph')
#   params:
#     win_size = 100000,
#     win_stride = 10000
#   threads: 12
#   conda: "../envs/hic_processing.yaml"
#   shell:
#     """
# 	tinycov covplot \
#       -N \
#       -n {wildcards.hic_library} \
#       -s {params.win_stride} \
#       -r {params.win_size} \
# 	    -t {output.text} \
#       -o {output.plot} \
#       {input}
#     """


# rule serpentine_binning:
#   input:
#     a = join(OUT, 'cool', "AT421.mcool"),
#     b = join(OUT, 'cool', "AT420.mcool")
#   output: join(OUT, 'plots', 'serpentine_i_u_ratio.svg')
#   params:
#     serp_res = LOW_RES
#   conda: '../envs/hic_processing.yaml'
#   shell:
#     """
#     python scripts/02_serpentine_analysis.py \
#       {input.a}::/resolutions/{params.serp_res} \
#       {input.b}::/resolutions/{params.serp_res} \
#       {output}
#     """

# rule compute_genomic_distance_law:
#   input: 
#     pairs = join(TMP, 'pairs', '{hic_library}.pairs'),
#     frags = join(TMP, 'digest', 'fragments_list.txt')
#   output:
#     tbl = join(TMP, 'distance_law', '{hic_library}_ps.tsv'),
#     plt = temp(join(TMP, 'distance_law', '{hic_library}_ps.svg'))
#   conda: '../envs/hic_processing.yaml'
#   shell:
#     """
#     hicstuff distancelaw -a \
#                          -i 1000 \
#                          -f {input.frags} \
#                          --pairs {input.pairs} \
#                          -O {output.tbl} \
#                          -o {output.plt}
#     """


# rule plot_distance_law:
#   input: expand(join(TMP, 'distance_law', '{hic_library}_ps.tsv'), hic_library=samples.library)
#   output: join(OUT, 'plots', 'distance_law_infection.svg')
#   params:
#     names = ','.join(samples.hic_library.values.tolist()),
#     inputs = ','.join([join(TMP, 'distance_law', f'{lib}_ps.tsv') for lib in samples.hic_library])
#   conda: '../envs/hic_processing.yaml'
#   shell: "hicstuff distancelaw -a -l {params.names} -o {output} --dist-tbl='{params.inputs}'"


# # Compute total contacts by condition for subsampling
# rule get_merged_contacts:
#   input:
#     lambda w: expand(
#       join(OUT, 'cool', '{hic_library}.cool'),
#       library=samples.loc[samples.condition==w.condition, 'hic_library']
#     )
#   output: join(TMP, '{condition}_contacts.txt')
#   conda: '../envs/hic_processing.yaml'
#   shell:
#     """
#     contacts=0
#     for cl in {input}
#     do
#       curr=$(cooler info $cl \
#         | grep sum \
#         | sed 's/[^0-9]*\([0-9]\+\)$/\1/')
#       contacts=$((contacts+curr))
#     done
#     echo $contacts > {output}
#     """


# # Retrieve the smallest number of contacts between conditions for subsampling
# rule get_lowest_contacts:
#   input: expand(join(TMP, '{condition}_contacts.txt'), condition=['infected', 'uninfected'])
#   output: join(TMP, 'lowest_condition_contacts.txt')
#   shell: "cat {input} | sort -k1,1n | head -1 > {output}"


# # Generate condition-merged cool files, subsampled at the same coverage for both conditions.
# rule condition_merged_cools:
#   input:
#     cool = lambda w: expand(
#       join(OUT, 'cool', '{hic_library}.cool'),
#       library=samples.loc[samples.condition==w.condition, 'library']
#     ),
#     target_contacts = join(TMP, 'lowest_condition_contacts.txt')
#   output: join(OUT, 'cool', 'sub_{condition}.mcool')
#   threads: 3
#   params:
#       out_cool = lambda w: join(OUT, 'cool', f'sub_{w.condition}.cool'),
#       max_res = MAX_RES,
#       med_res = MED_RES,
#       low_res = LOW_RES
#   conda: "../envs/hic_processing.yaml"
#   shell:
#     """
#     COUNT=$(head -n1 {input.target_contacts})
#     cooler merge {params.out_cool}_tmp {input.cool}
#     # Skips subsampling for the condition with lowest coverage
#     cooltools random-sample --count "$COUNT" {params.out_cool}_tmp {params.out_cool} \
#       || cp {params.out_cool}_tmp {params.out_cool}
# 		cooler zoomify -r {params.max_res},{params.med_res},{params.low_res} \
# 		  --balance \
# 		  --balance-args '--mad-max 10' \
# 		  -n {threads} \
# 		  {params.out_cool}
#     rm {params.out_cool} {params.out_cool}_tmp
#     """