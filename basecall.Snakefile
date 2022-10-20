
def fix_name(new_name):
    """
    Terrible hack. Sets the name of the most recently created rule to be
    new_name.
    """
    list(workflow.rules)[-1].name = new_name
    temp_rules = list(rules.__dict__.items())
    temp_rules[-1] = (new_name, temp_rules[-1][1]) 
    rules.__dict__ = dict(temp_rules)


# Local fast5 files
fast5_path = 'data/reads/BB31_drone'
# fast5_path = 'data/reads/basecalling_practical' # from https://timkahlke.github.io/LongRead_tutorials

versions_to_run = [
    'guppy_3.4.4',
    'guppy_3.6.0',
    'guppy_4.0.14',
    'guppy_4.2.2',
    'guppy_4.5.4',
    'guppy_5.0.16',
    'guppy_5.0.16_sup',
    'guppy_6.1.3',
    'guppy_6.1.3_sup']
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
            partition = 'gpgpu',
            qos = 'gpgpumdhs',
            gres = 'gpu:1',
            proj = 'punim1712',
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
