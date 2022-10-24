#!/usr/bin/env python3

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from pathlib import Path
import os


#############
# FUNCTIONS #
#############

def busco_input(wildcards):
    if wildcards.guppy == 'ref':
        return(raw_ref)
    else:
        return('output/051_oriented/{guppy}.{flye_mode}/contigs.fa')


def combine_indiv_reads(wildcards):
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
    my_output_path = f'output/tmp/020_porechop/{wildcards.guppy}/{{read}}.fastq'
    my_read_names = snakemake.io.glob_wildcards(my_read_path).read
    # check file size
    non_empty_read_names = (
        x for x in my_read_names if os.stat(my_read_path.format(read=x)).st_size > 0)
    my_output = snakemake.io.expand(my_output_path, read=non_empty_read_names)
    return(sorted(set(my_output)))
    # return(sorted(set(x for x in my_reads if os.stat(x).st_size > 0)))
    

def fix_name(new_name):
    """
    Terrible hack. Sets the name of the most recently created rule to be
    new_name.
    """
    list(workflow.rules)[-1].name = new_name
    temp_rules = list(rules.__dict__.items())
    temp_rules[-1] = (new_name, temp_rules[-1][1]) 
    rules.__dict__ = dict(temp_rules)


def get_porechop_input(wildcards):
    # need to handle old version of guppy that don't have pass/fail dirs
    if Path(f'output/010_basecall/{wildcards.guppy}/pass').is_dir():
        return(f'output/010_basecall/{{guppy}}/pass/{{read}}.fastq')
    else:
        return(f'output/010_basecall/{{guppy}}/{{read}}.fastq')


###########
# GLOBALS #
###########

# NCBI reference genome
HTTP = HTTPRemoteProvider()
remote_ref_url = (
    'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/'
    'GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz')

local_ref = 'output/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz'
raw_ref = 'output/GCF_003254395.2_Amel_HAv3.1_genomic.fna'


# remote_ref = HTTP.remote(
#     ('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/405/'
#      'GCF_000149405.2_ASM14940v2/GCF_000149405.2_ASM14940v2_assembly_structure/'
#      'Primary_Assembly/assembled_chromosomes/FASTA/chr17.fna.gz'),
#     keep_local=True)
# raw_ref = 'data/GCF_000149405.2_chr17.fna'

# BUSCO lineage
busco_lineage = 'hymenoptera_odb10'
lineage_url = ('https://busco-data.ezlab.org/v5/data/lineages/'
               'hymenoptera_odb10.2020-08-05.tar.gz')
lineage_archive = f'output/080_busco/{busco_lineage}.tar.gz'

# busco_lineage = 'eukaryota_odb10'
# lineage_archive = HTTP.remote(
#     ('https://busco-data.ezlab.org/v5/data/lineages/'
#      'eukaryota_odb10.2020-09-10.tar.gz'),
#     keep_local=True)
lineage_path = f'output/080_busco/{busco_lineage}'


# guppy version I have
# sync with basecall snakefile
versions_to_run = [
    'guppy_3.4.4',
    'guppy_3.6.0',
    'guppy_4.0.14',
    'guppy_4.2.2',
    'guppy_4.5.4',
    'guppy_5.0.16',
    'guppy_5.0.16_sup',
    'guppy_6.1.3',
    'guppy_6.1.3_sup',
    'guppy_6.3.8',
    'guppy_6.3.8_sup']
guppy_versions = {
    'guppy_3.4.1': 'shub://TomHarrop/ont-containers:guppy_3.4.1',
    'guppy_3.4.4': 'shub://TomHarrop/ont-containers:guppy_3.4.4',
    'guppy_3.6.0': 'shub://TomHarrop/ont-containers:guppy_3.6.0',
    'guppy_4.0.11': 'shub://TomHarrop/ont-containers:guppy_4.0.11',
    'guppy_4.0.14': 'shub://TomHarrop/ont-containers:guppy_4.0.14',
    'guppy_4.2.2': 'shub://TomHarrop/ont-containers:guppy_4.2.2',
    'guppy_4.5.4': 'docker://ghcr.io/tomharrop/container-guppy:4.5.4',
    'guppy_5.0.16': 'docker://ghcr.io/tomharrop/container-guppy:5.0.16',
    'guppy_5.0.16_sup': 'docker://ghcr.io/tomharrop/container-guppy:5.0.16',
    'guppy_6.1.3': 'docker://ghcr.io/tomharrop/container-guppy:6.1.3',
    'guppy_6.1.3_sup': 'docker://ghcr.io/tomharrop/container-guppy:6.1.3', # dna_r9.4.1_450bps_sup.cfg
    'guppy_6.3.8': 'docker://ghcr.io/tomharrop/container-guppy:6.3.8',
    'guppy_6.3.8_sup': 'docker://ghcr.io/tomharrop/container-guppy:6.3.8' # dna_r9.4.1_450bps_sup.cfg
}

# Containers
biopython = 'docker://quay.io/biocontainers/biopython:1.78'
bcftools = 'docker://quay.io/biocontainers/bcftools:1.15.1--h0ea216a_0'
busco = 'docker://quay.io/biocontainers/busco:5.3.2--pyhdfd78af_0'
filtlong = 'docker://quay.io/biocontainers/filtlong:0.2.1--hd03093a_1'
flye = 'docker://quay.io/biocontainers/flye:2.9--py39h6935b12_1'
minimap = 'docker://quay.io/biocontainers/minimap2:2.24--h7132678_1'
mummer = 'docker://quay.io/biocontainers/mummer:3.23--pl5321h1b792b2_13'
porechop = 'docker://quay.io/biocontainers/porechop:0.2.4--py39hc16433a_3'
quast = 'docker://quay.io/biocontainers/quast:5.0.2--py36pl5321hcac48a8_7'
ragtag = 'docker://quay.io/biocontainers/ragtag:2.1.0--pyhb7b1952_0'
samtools = 'docker://quay.io/biocontainers/samtools:1.15.1--h1170115_0'


#########
# RULES #
#########


subworkflow basecall:
    snakefile: 'basecall.Snakefile'

wildcard_constraints:
    guppy = '|'.join(versions_to_run) + '|ref'

rule target:
    input:
        # expand('output/060_dnadiff/{guppy}.{flye_mode}/contigs.snps',
        #        guppy=versions_to_run,
        #        flye_mode=['nano-raw', 'nano-hq']),
        expand('output/065_minimap-snps/{guppy}.{flye_mode}/out.snps.vcf',
               guppy=versions_to_run,
               flye_mode=['nano-raw', 'nano-hq']),
        expand('output/080_busco/{guppy}.{flye_mode}/run_{busco_lineage}/full_table.tsv',
               busco_lineage=busco_lineage,
               guppy=versions_to_run,
               flye_mode=['nano-raw', 'nano-hq']),
        f'output/080_busco/ref.ref/run_{busco_lineage}/full_table.tsv',
        'output/tmp/070_quast/report.txt'


# compare genomes
# using busco
rule busco:
    input:
        fasta = busco_input,
        lineage = lineage_path
    output:
        f'output/080_busco/{{guppy}}.{{flye_mode}}/run_{busco_lineage}/full_table.tsv',
    log:
        Path(('output/logs/'
              'busco.{guppy}.{flye_mode}.log')).resolve()
    benchmark:
        Path('output/benchmarks/busco.{guppy}.{flye_mode}.log.tsv').resolve()
    params:
        wd = 'output/080_busco',
        fasta = lambda wildcards, input:
            Path(input.fasta).resolve(),
        lineage = lambda wildcards, input:
            Path(input.lineage).resolve()
    threads:
        workflow.cores
    resources:
        time = 60,
        mem_mb = 10000
    singularity:
        busco
    shell:
        'mkdir -p {params.wd} ; '
        'cd {params.wd} || exit 1 ; '
        'busco '
        '--offline '
        '--force '
        '--in {params.fasta} '
        '--out {wildcards.guppy}.{wildcards.flye_mode} '
        '--lineage_dataset {params.lineage} '
        '--cpu {threads} '
        '--mode genome '
        '&> {log}'

# using quast
rule quast:
    input:
        genomes = expand('output/051_oriented/{guppy}.{flye_mode}/contigs.fa',
                         guppy=versions_to_run,
                         flye_mode=['nano-raw', 'nano-hq']),
        ref = local_ref
    output:
        'output/tmp/070_quast/report.txt'
    params:
        outdir = 'output/tmp/070_quast'
    threads:
        workflow.cores
    resources:
        time = 1440,
        mem_mb = 50000
    log:
        'output/logs/quast.log'
    benchmark:
        'output/benchmarks/quast.tsv'
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


# needs to be replaced with minimap2
# docker://quay.io/biocontainers/minimap2:2.24--h7132678_1 has paftools
# see https://github.com/lh3/minimap2/blob/fe35e679e95d936698e9e937acc48983f16253d6/cookbook.md#calling-variants-from-assembly-to-reference-alignment
rule snps_only:
    input:
        'output/065_minimap-snps/{guppy}.{flye_mode}/out.vcf'
    output:
        'output/065_minimap-snps/{guppy}.{flye_mode}/out.snps.vcf'
    log:
        'output/logs/snps_only.{guppy}.{flye_mode}.log'
    container:
        bcftools
    shell:
        'bcftools view '
        '--types snps '
        '{input} '
        '> {output} '
        '2> {log}'


rule paftools_snps:
    input:
        ref = raw_ref,
        paf = 'output/tmp/065_minimap-snps/{guppy}.{flye_mode}/aln.sorted.paf'
    output:
        'output/065_minimap-snps/{guppy}.{flye_mode}/out.vcf'
    log:
        paftools = 'output/logs/minimap_snps.{guppy}.{flye_mode}.paftools.log',
    resources:
        time = 10,
    container:
        minimap
    shell:
        'paftools.js '
        'call '
        '-f {input.ref} '
        '<( cat {input.paf} ) '
        '>> {output} '
        '2> {log}'


rule minimap_sort:
    input:
        'output/tmp/065_minimap-snps/{guppy}.{flye_mode}/aln.paf'
    output:
        pipe('output/tmp/065_minimap-snps/{guppy}.{flye_mode}/aln.sorted.paf')
    threads:
        1
    resources:
        time = 10
    container:
        minimap
    shell:
        'sort -k6,6 -k8,8n <( cat {input} ) >> {output} '

rule minimap:
    input:
        ref = raw_ref,
        contigs = 'output/051_oriented/{guppy}.{flye_mode}/contigs.fa',
    output:
        pipe('output/tmp/065_minimap-snps/{guppy}.{flye_mode}/aln.paf')
    log:
        'output/logs/minimap_snps.{guppy}.{flye_mode}.minimap.log',
    benchmark:
        'output/benchmarks/minimap_snps.{guppy}.{flye_mode}.minimap.tsv'
    threads:
        min(12, workflow.cores - 2)
    resources:
        time = 10,
        mem_mb = 4000 * 4
    container:
        minimap
    shell:
        'minimap2 '
        '-t  {threads} '
        '-cx asm5 '
        '--cs '
        '{input.ref} '
        '{input.contigs} '
        '>> {output} '
        '2> {log} '


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
    resources:
        time = 2880,
        mem_mb = 50000
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
        ref = local_ref,
        query = 'output/040_flye/{guppy}.{flye_mode}/assembly.fasta'
    output:
        'output/050_ragtag/{guppy}.{flye_mode}/ragtag.scaffold.fasta',
        'output/050_ragtag/{guppy}.{flye_mode}/ragtag.scaffold.agp'
    params:
        wd = 'output/050_ragtag/{guppy}.{flye_mode}'
    log:
        'output/logs/ragtag.{guppy}.{flye_mode}.log'
    benchmark:
        'output/benchmarks/ragtag.{guppy}.{flye_mode}.tsv'
    threads:
        min(workflow.cores, 64)
    resources:
        time = 120,
        mem_mb = 50000
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
        fq = 'output/tmp/030_filtlong/{guppy}.fastq'
    output:
        'output/040_flye/{guppy}.{flye_mode}/assembly.fasta'
    params:
        outdir = 'output/040_flye/{guppy}.{flye_mode}',
        mode = '{flye_mode}'
    threads:
        min(128, workflow.cores)
    resources:
        time = 120,
        mem_mb = 128000
    log:
        'output/logs/flye.{guppy}.{flye_mode}.log'
    benchmark:
        'output/benchmarks/flye.{guppy}.{flye_mode}.tsv'
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
        'output/tmp/020_porechop/{guppy}.fastq'
    output:
        'output/tmp/030_filtlong/{guppy}.fastq'
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

rule combine_indiv_reads:
    input:
        combine_indiv_reads
    output:
        'output/tmp/020_porechop/{guppy}.fastq'
    resources:
        time = 10
    container:
        porechop
    shell:
        'cat {input} > {output}'


rule porechop:
    input:
        basecall('output/010_basecall/{guppy}/sequencing_summary.txt'),
        read = get_porechop_input
    output:
        'output/tmp/020_porechop/{guppy}/{read}.fastq'
    log:
        'output/logs/porechop.{guppy}.{read}.log'
    threads:
        1
    resources:
        time = 10
    container:
        porechop
    shell:
        'porechop '
        '-i {input.read} '
        '-o {output} '
        '--verbosity 1 '
        '--threads {threads} '
        '--discard_middle '
        '&> {log}'


# GENERIC
rule busco_expand:
    input:
        lineage_archive
    output:
        directory(lineage_path)
    singularity:
        busco
    shell:
        'mkdir -p {output} && '
        'tar -zxf '
        '{input} '
        '-C {output} '
        '--strip-components 1 '

rule download_lineage:
    input:
        HTTP.remote(
            lineage_url,
            keep_local=True)
    output:
        lineage_archive
    shell:
        'mv {input} {output}'

rule raw_ref:
    input:
        local_ref
    output:
        fa = raw_ref,
        fai = f'{raw_ref}.fai'
    singularity:
        samtools
    shell:
        'gunzip -c {input} > {output.fa} ; '
        'samtools faidx {output.fa}'

rule download_ref:
    input:
        HTTP.remote(
            remote_ref_url,
            keep_local=True)
    output:
        local_ref
    shell:
        'mv {input} {output}'
