#!/usr/bin/env python3

#############
# FUNCTIONS #
#############

def choose_container(guppy_version):
    return(guppy_versions[guppy_version])


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


# drop reads < 5kb
# remove worst 10% of reads (check cov)
# get IDs
# only basecall those (use option -l in guppy)



#########
# RULES #
#########


rule target:
    input:
        'output/010_basecall/full/guppy_6.1.3/sequencing_summary.txt'


# rule porechop:
#     input:
#         'output/010_basecall/full/guppy_6.1.3/sequencing_summary.txt'
#     output:
#         temp('output/010_porechop/{indiv}/{read}.fastq')
#     log:
#         'output/logs/porechop.{indiv}.{read}.log'
#     threads:
#         1
#     container:
#         porechop
#     shell:
#         'porechop '
#         '-i {input} '
#         '-o {output} '
#         '--verbosity 1 '
#         '--threads {threads} '
#         '--discard_middle '
#         '&> {log}'



rule full_basecall:
    input:
        'data/reads/BB31_drone'
    output:
        'output/010_basecall/full/guppy_6.1.3/sequencing_summary.txt'
    params:
        outdir = 'output/010_basecall/full/guppy_6.1.3',
        flowcell = "FLO-MIN106",
        kit = "SQK-LSK109"
    log:
        'output/logs/full_basecall.guppy_6.1.3.log'
    container:
        guppy_versions['guppy_6.1.3']
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


