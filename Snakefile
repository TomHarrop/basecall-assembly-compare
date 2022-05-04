#!/usr/bin/env python3

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider


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

# NCBI reference genome
HTTP = HTTPRemoteProvider()

# Amel_HAv3 = HTTP.remote(
#     ('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/'
#      'GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz'),
#     keep_local=True)
# raw_ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'

remote_ref = HTTP.remote(
    ('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/405/'
     'GCF_000149405.2_ASM14940v2/GCF_000149405.2_ASM14940v2_assembly_structure/'
     'Primary_Assembly/assembled_chromosomes/FASTA/chr17.fna.gz'),
    keep_local=True)
raw_ref = 'data/GCF_000149405.2_chr17.fna'

# Local fast5 files
# fast5_path = 'data/reads/BB31_drone'
fast5_path = 'data/reads/basecalling_practical' # from https://timkahlke.github.io/LongRead_tutorials

# guppy version I have
versions_to_run = ['guppy_6.1.3', 'guppy_3.6.0']
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

# Containers
biopython = 'docker://quay.io/biocontainers/biopython:1.78'
filtlong = 'docker://quay.io/biocontainers/filtlong:0.2.1--hd03093a_1'
flye = 'docker://quay.io/biocontainers/flye:2.9--py39h6935b12_1'
mummer = 'docker://quay.io/biocontainers/mummer:3.23--pl5321h1b792b2_13'
porechop = 'docker://quay.io/biocontainers/porechop:0.2.4--py39hc16433a_3'
ragtag = 'docker://quay.io/biocontainers/ragtag:2.1.0--pyhb7b1952_0'
quast = 'docker://quay.io/biocontainers/quast:5.0.2--py36pl5321hcac48a8_7'


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
        expand('output/060_dnadiff/{guppy}.{flye_mode}/contigs.snps',
               guppy=versions_to_run,
               flye_mode=['nano-raw', 'nano-hq']),
        'output/070_quast/report.txt'





# compare genomes
# using quast
rule quast:
    input:
        genomes = expand('output/051_oriented/{guppy}.{flye_mode}/contigs.fa',
                         guppy=versions_to_run,
                         flye_mode=['nano-raw', 'nano-hq']),
        ref = remote_ref
    output:
        'output/070_quast/report.txt'
    params:
        outdir = 'output/070_quast'
    threads:
        workflow.cores
    log:
        'output/logs/quast.log'
    container:
        quast
    shell:
        'quast '
        '-o {params.outdir} '
        '-r {input.ref} '
        '-t {threads} '
        '-L '
        '--eukaryote '
        '--k-mer-stats '
        '{input} '
        '&> {log}'


# dnadiff uses Olin Silander's method (https://github.com/osilander/bonito_benchmarks)
rule dnadiff:
    input:
        ref = raw_ref,
        contigs = 'output/051_oriented/{guppy}.{flye_mode}/contigs.fa',
    output:
        'output/060_dnadiff/{guppy}.{flye_mode}/contigs.snps'
    params:
        prefix = "output/060_dnadiff/{guppy}.{flye_mode}/contigs",
    log:
        'output/logs/dnadiff.{guppy}.{flye_mode}.log'
    container:
        mummer
    shell:
        'dnadiff '
        '{input.ref} '
        '{input.contigs} '
        '-p {params.prefix} '
        '&> {log}'


# extract genomic contigs.
# n.b. this removes unplaced contigs!!!
rule orient_scaffolds:
    input:
        fa = 'output/040_flye/{guppy}.{flye_mode}/assembly.fasta',
        agp = 'output/050_ragtag/{guppy}.{flye_mode}/ragtag.scaffold.agp'
    output:
        fa = 'output/051_oriented/{guppy}.{flye_mode}/contigs.fa'
    log:
        'output/logs/orient_scaffolds.{guppy}.{flye_mode}.log'
    container:
        biopython
    script:
        'src/orient_scaffolds.py'

rule ragtag:
    input:
        ref = remote_ref,
        query = 'output/040_flye/{guppy}.{flye_mode}/assembly.fasta'
    output:
        'output/050_ragtag/{guppy}.{flye_mode}/ragtag.scaffold.fasta',
        'output/050_ragtag/{guppy}.{flye_mode}/ragtag.scaffold.agp'
    params:
        wd = 'output/050_ragtag/{guppy}.{flye_mode}'
    log:
        'output/logs/ragtag.{guppy}.{flye_mode}.log'
    threads:
        min(workflow.cores, 64)
    container:
        ragtag
    shell:
        'ragtag.py scaffold '
        '-o {params.wd} '
        '-w '
        # '-r -g 101 '    # only add gaps 101 Ns or longer DOESN'T WORK
        '-t {threads} '
        '{input.ref} '
        '{input.query} '
        '&> {log}'



# run assembly
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
        '--target_bases 50000000 ' # this is almost 100x for diatom
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
            fast5_path
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

    # fix_name(guppy)

# GENERIC
rule raw_ref:
    input:
        remote_ref
    output:
        raw_ref
    singularity:
        flye
    shell:
        'gunzip -c {input} > {output}'