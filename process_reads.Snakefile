#!/usr/bin/env python3

from pathlib import Path
import csv

#############
# FUNCTIONS #
#############

def get_guppy_fastq_files(wildcards):
    # need to handle old version of guppy that don't have pass/fail dirs
    if Path(f'output/010_basecall/{wildcards.guppy}/pass').is_dir():
        return(f'output/010_basecall/{{guppy}}/pass/{{read}}.fastq')
    else:
        return(f'output/010_basecall/{{guppy}}/{{read}}.fastq')


def aggregate_reads(wildcards):
    # Hack: mark ALL the basecall checkpoints as required for porechop. This
    # prevents the glob from happening until all the basecall jobs are
    # finished
    # for it in vars(checkpoints):
    #     try:
    #         checkpoints.__dict__[it].get()
    #     except AttributeError:
    #         pass
    # bail if the summary file isn't there
    summary_file = f'output/010_basecall/{wildcards.guppy}/sequencing_summary.txt'
    try:
        os.stat(summary_file)
    except FileNotFoundError:
        print(f'ERROR. {summary_file} not found.')
        print('       Run basecall.Snakefile')
        raise FileNotFoundError
    # need to handle old version of guppy that don't have pass/fail dirs
    if Path(f'output/010_basecall/{wildcards.guppy}/pass').is_dir():
        my_read_path = f'output/010_basecall/{wildcards.guppy}/pass/{{read}}.fastq'
    else:
        my_read_path = f'output/010_basecall/{wildcards.guppy}/{{read}}.fastq'
    my_output_path = f'output/010_basecall/{wildcards.guppy}/compressed_reads/{{read}}.fastq.gz'
    my_read_names = snakemake.io.glob_wildcards(my_read_path).read
    # check file size
    non_empty_read_names = (
        x for x in my_read_names if os.stat(my_read_path.format(read=x)).st_size > 0)
    my_output = snakemake.io.expand(my_output_path, read=non_empty_read_names)
    return(sorted(set(my_output)))
    # return(sorted(set(x for x in my_reads if os.stat(x).st_size > 0)))


###########
# GLOBALS #
###########

versions_manifest = 'data/versions_to_run.csv'

# CONTAINERS
pigz = 'docker://quay.io/biocontainers/pigz:2.3.4'

########
# MAIN #
########

guppy_versions = {}
reader = csv.DictReader(open(versions_manifest, 'rt'))
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
		expand('output/test/reads/{guppy}.fastq.gz',
			   guppy=['guppy_4.5.4'])

rule aggregate_reads:
	input:
		aggregate_reads
	output:
		'output/test/reads/{guppy}.fastq.gz'
	shell:
		'cat {input} > {output}'

# compress the guppy output
rule gzip_fastq_files:
    input:
        basecall('output/010_basecall/{guppy}/sequencing_summary.txt'),
        # this has wildcards {guppy} and {read}
        read = get_guppy_fastq_files
    output:
        'output/010_basecall/{guppy}/compressed_reads/{read}.fastq.gz'
    log:
        'output/logs/gzip_fastq_files/{guppy}.{read}.log'
    threads:
        10
    resources:
        time = 1
    container:
    	pigz
    shell:
        'pigz -9 <{input.read} >{output} '
        '&& rm {input.read} '
        '&> {log}'
