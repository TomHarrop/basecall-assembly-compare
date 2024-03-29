#!/usr/bin/env python3

from pathlib import Path
import csv

#############
# FUNCTIONS #
#############

def get_guppy_fastq_files(wildcards):
    #  make sure basecalling has happened
    summary_file = f'output/010_basecall/{wildcards.guppy}/sequencing_summary.txt'
    try:
        os.stat(summary_file)
    except FileNotFoundError:
        print(f'ERROR. {summary_file} not found.')
        print('       Run basecall.Snakefile')
        raise FileNotFoundError
    # need to handle old version of guppy that don't have pass/fail dirs
    if Path(f'output/010_basecall/{wildcards.guppy}/pass').is_dir():
        return(f'output/010_basecall/{{guppy}}/pass/{{read}}.fastq')
    else:
        return(f'output/010_basecall/{{guppy}}/{{read}}.fastq')


def aggregate_reads(wildcards):
    idlist = checkpoints.generate_read_id_list.get(**wildcards).output['idlist']
    with open(idlist, 'rt') as f:
        read_ids = [line.rstrip() for line in f]
    return(
        snakemake.io.expand(
            'output/tmp/020_porechop/{{guppy}}/{read}.fastq',
            read=read_ids))

###########
# GLOBALS #
###########

versions_manifest = 'data/versions_to_run.csv'

# CONTAINERS
biopython = 'docker://quay.io/biocontainers/biopython:1.78'
filtlong = 'docker://quay.io/biocontainers/filtlong:0.2.1--hd03093a_1'
pigz = 'docker://quay.io/biocontainers/pigz:2.3.4'
porechop = 'docker://quay.io/biocontainers/porechop:0.2.4--py39hc16433a_3'

########
# MAIN #
########

guppy_versions = {}
with open(versions_manifest, 'rt') as f:
    reader = csv.DictReader(row for row in f if not row.startswith('#'))
    for row in reader:
        guppy_versions[row['name']] = row['container']
    
versions_to_run = sorted(set(guppy_versions.keys()))


#########
# RULES #
#########

subworkflow basecall:
    snakefile: 'basecall.Snakefile'

wildcard_constraints:
    guppy = '|'.join(versions_to_run) + '|ref'

rule target:
    input:
        expand('output/035_processed_reads/{guppy}.fastq.gz',
               guppy=versions_to_run)

rule compress_processed_reads:
    input:
        'output/tmp/030_filtlong/{guppy}.fastq'
    output:
        'output/035_processed_reads/{guppy}.fastq.gz'
    log:
        'output/logs/compress_processed_reads/{guppy}.log'
    threads:
        10
    resources:
        time = 20
    container:
        pigz
    shell:
        'pigz -p {threads} '
        '-9 '
        '<{input} '
        '>{output} '
        '2> {log}'

rule filtlong:
    input:
        'output/tmp/020_porechop/{guppy}.fastq'
    output:
        temp('output/tmp/030_filtlong/{guppy}.fastq')
    log:
        'output/logs/filtlong.{guppy}.log'
    resources:
        time = 120
    container:
        filtlong
    shell:
        'filtlong '
        # '--target_bases 50000000 ' # this is almost 100x for diatom
        '--target_bases 8000000000 ' # 10 GB is approx 50x for amel
        '--min_length 5000 '
        '{input} '
        '> {output} '
        '2> {log}'

rule aggregate_reads:
    input:
        aggregate_reads
    output:
        temp('output/tmp/020_porechop/{guppy}.fastq')
    shell:
        'cat {input} > {output}'

rule porechop:
    input:
        'output/tmp/010_basecall/{guppy}/{read}.fastq'
    output:
        temp('output/tmp/020_porechop/{guppy}/{read}.fastq')
    log:
        'output/logs/porechop/{guppy}.{read}.log'
    threads:
        1
    resources:
        time = 10
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


# extract the guppy output for processing
rule unzip_fastq_file:
    input:
        'output/010_basecall/{guppy}/compressed_reads/{read}.fastq.gz'
    output:
        temp('output/tmp/010_basecall/{guppy}/{read}.fastq')
    log:
        'output/logs/unzip_fastq_file/{guppy}.{read}.log'
    container:
        pigz
    shell:
        'pigz -d <{input} >{output} 2>{log}'

# compress the guppy output for storage
rule gzip_fastq_file:
    input:
        read = get_guppy_fastq_files
    output:
        'output/010_basecall/{guppy}/compressed_reads/{read}.fastq.gz'
    log:
        'output/logs/gzip_fastq_file/{guppy}.{read}.log'
    threads:
        10
    resources:
        time = 1
    container:
        pigz
    shell:
        'pigz -p {threads} '
        '-9 '
        '<{input.read} '
        '>{output} '
        '&& rm {input.read} '
        '&> {log}'

# generate a list of read IDs
checkpoint generate_read_id_list:
    input:
        seqsum = basecall('output/010_basecall/{guppy}/sequencing_summary.txt')
    output:
        idlist = 'output/011_read_ids/{guppy}/read_list.txt'
    threads:
        1
    resources:
        time = 1
    log:
        'output/logs/generate_read_id_list/{guppy}.log'
    container:
        biopython
    script:
        'src/generate_read_id_list.py'
