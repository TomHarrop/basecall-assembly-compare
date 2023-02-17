#!/usr/bin/env python3

import csv

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

versions_manifest = 'data/versions_to_run.csv'
# Local fast5 files
fast5_path = 'data/reads/BB31_drone'
# fast5_path = 'data/reads/basecalling_practical' # from https://timkahlke.github.io/LongRead_tutorials


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

wildcard_constraints:
    guppy = '|'.join(versions_to_run)


rule target:
    input:
        expand('output/010_basecall/{guppy}/sequencing_summary.txt',
               guppy=versions_to_run)



# Full basecall. Could reduce with this strategy:
# drop reads < 5kb
# remove worst 10% of reads (check cov)
# get IDs
# only basecall those (use option -l in guppy)
for guppy in versions_to_run:
    rule:
        input:
            fast5_path
        output:
            f'output/010_basecall/{guppy}/sequencing_summary.txt',
            # p = directory(f'output/010_basecall/{guppy}/pass'),
            # f = directory(f'output/010_basecall/{guppy}/fail')
        params:
            outdir = f'output/010_basecall/{guppy}',
            config = ('dna_r9.4.1_450bps_sup.cfg' if guppy.endswith('_sup')
                      else 'dna_r9.4.1_450bps_hac.cfg')
        log:
            f'output/logs/full_basecall.{guppy}.log'
        threads:
            3
        resources:
            partition = 'gpu-a100',
            gres = 'gpu:1',
            time = 480 * 5,
            mem_mb = 40000
        container:
            guppy_versions[guppy]
        shell:
            # 'nvidia-smi && '
            'guppy_basecaller '
            '--device auto '        # enable GPU
            '--input_path {input} '
            '--save_path {params.outdir} '
            '--config {params.config} '
            '--verbose_logs '
            '--recursive '
            '&> {log}'

    fix_name(guppy)
