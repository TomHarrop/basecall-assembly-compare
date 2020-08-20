#!/usr/bin/env python3


#############
# FUNCTIONS #
#############

def combine_indiv_reads(wildcards):
    # Hack: mark ALL the basecall checkpoints as required for porechop. This
    # prevents the glob from happening until all the basecall jobs are
    # finished
    for it in vars(checkpoints):
        try:
            checkpoints.__dict__[it].get()
        except AttributeError:
            pass
    my_read_path = f'output/010_basecall/{wildcards.guppy}/pass/{{read}}.fastq'
    my_output_path = f'output/020_porechop/{wildcards.guppy}/{{read}}.fastq'
    my_read_names = snakemake.io.glob_wildcards(my_read_path).read
    my_output = snakemake.io.expand(my_output_path, read=my_read_names)
    return(sorted(set(my_output)))


def fix_name(new_name):
    """
    Terrible hack. Sets the name of the most recently created rule to be
    new_name.
    """
    list(workflow.rules)[-1].name = new_name
    temp_rules = list(rules.__dict__.items())
    temp_rules[-1] = (new_name, temp_rules[-1][1]) 
    rules.__dict__ = dict(temp_rules)


###########
# GLOBALS #
###########

versions_to_run = ['guppy_6.1.3']

guppy_versions = {
    'guppy_3.4.1': 'shub://TomHarrop/ont-containers:guppy_3.4.1',
    'guppy_3.4.4': 'shub://TomHarrop/ont-containers:guppy_3.4.4',
    'guppy_3.6.0': 'shub://TomHarrop/ont-containers:guppy_3.6.0',
    'guppy_4.0.11': 'shub://TomHarrop/ont-containers:guppy_4.0.11',
    'guppy_4.0.14': 'shub://TomHarrop/ont-containers:guppy_4.0.14',
    'guppy_4.2.2': 'shub://TomHarrop/ont-containers:guppy_4.2.2',
    'guppy_4.5.4': 'docker://ghcr.io/tomharrop/container-guppy:4.5.4',
    'guppy_5.0.16': 'docker://ghcr.io/tomharrop/container-guppy:5.0.16',
    'guppy_6.1.3': 'docker://ghcr.io/tomharrop/container-guppy:6.1.3',
}


filtlong = 'docker://quay.io/biocontainers/filtlong:0.2.1--hd03093a_1'
flye = 'docker://quay.io/biocontainers/flye:2.9--py39h6935b12_1'
porechop = 'docker://quay.io/biocontainers/porechop:0.2.4--py39hc16433a_3'

# drop reads < 5kb
# remove worst 10% of reads (check cov)
# get IDs
# only basecall those (use option -l in guppy)



#########
# RULES #
#########

wildcard_constraints:
    guppy = '|'.join(versions_to_run),

rule target:
    input:
        expand('output/040_flye/{guppy}.{flye_mode}/assembly.fasta',
               guppy=versions_to_run,
               flye_mode=['nano-raw', 'nano-hq']),


rule flye:
    input:
        fq = 'output/030_filtlong/{guppy}.fastq'
    output:
        'output/040_flye/{guppy}.{flye_mode}/assembly.fasta'
    params:
        outdir = 'output/040_flye/{guppy}.{flye_mode}',
        mode = '{flye_mode}'
    threads:
        min(128, workflow.cores)
    log:
        'output/logs/flye.{guppy}.{flye_mode}.log'
    container:
        flye
    shell:
        'flye '
        # '--resume '
        '--{wildcards.flye_mode} '
        '{input.fq} '
        '--out-dir {params.outdir} '
        '--threads {threads} '
        '&>> {log}'


# PROCESS READS
rule filtlong:
    input:
        'output/020_porechop/{guppy}.fastq'
    output:
        'output/030_filtlong/{guppy}.fastq'
    log:
        'output/logs/filtlong.{guppy}.log'
    container:
        filtlong
    shell:
        'filtlong '
        '--target_bases 1500000000 ' # this is approx 10x for testing the pipeline
        '--min_length 5000 '
        '{input} '
        '> {output} '
        '2> {log}'

rule combine_indiv_reads:
    input:
        combine_indiv_reads
    output:
        'output/020_porechop/{guppy}.fastq'
    container:
        porechop
    shell:
        'cat {input} > {output}'

rule porechop:
    input:
        'output/010_basecall/{guppy}/pass/{read}.fastq'
    output:
        temp('output/020_porechop/{guppy}/{read}.fastq')
    log:
        'output/logs/porechop.{guppy}.{read}.log'
    threads:
        1
    container:
        porechop
    shell:
        'porechop '
        '-i {input} '
        '-o {output} '
        '--verbosity 1 '
        '--threads {threads} '
        '--discard_middle '
        '&> {log}'

for guppy in versions_to_run:
    checkpoint:
        input:
            'data/reads/BB31_drone'
        output:
            f'output/010_basecall/{guppy}/sequencing_summary.txt',
            p = directory(f'output/010_basecall/{guppy}/pass'),
            f = directory(f'output/010_basecall/{guppy}/fail')
        params:
            outdir = f'output/010_basecall/{guppy}',
            flowcell = "FLO-MIN106",
            kit = "SQK-LSK109"
        log:
            f'output/logs/full_basecall.{guppy}.log'
        container:
            guppy_versions[guppy]
        shell:
            'guppy_basecaller '
            '--device auto '        # enable GPU
            '--input_path {input} '
            '--save_path {params.outdir} '
            '--flowcell {params.flowcell} '
            '--kit {params.kit} '
            '--verbose_logs '
            '--recursive '
            '&> {log}'

    fix_name(guppy)
