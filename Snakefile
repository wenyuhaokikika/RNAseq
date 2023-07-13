'''
@File    :   Snakefile
@Time    :   2023/01/30 14:23:11
@Author  :   wenyuhao 
@Version :   1.0
@Contact :   925201392@qq.com
@Desc    :   None
'''

# here put the import lib
include: 'rules/common.smk'
include: 'rules/ref.smk'
include: 'rules/preprocess.smk'
include: 'rules/mapping.smk'
include: 'rules/aggregate.smk'
#include: 'rules/comparison.smk'


import os

BASE_DIR = os.path.dirname(workflow.snakefile)


rule all:
    input:
        # rules.genome.output.fasta,
        # rules.genome.output.gtf,
        # rules.sequenceDownload.output.cdna,
        # rules.sequenceDownload.output.pep,
        # rules.sequenceDownload.output.cds,
        rules.STARindex.output,
        rules.annotation.output,
        [refPath(i) for i in ['genome.fa.gz','genome.gtf.gz','cdna.fa.gz','pep.fa.gz','cds.fa.gz']],
        # [directory(refPath(i)) for i in ['STARindex','annotation']],
        expand(
            config['workspace'] + '/aggregate/all_sample_{set}_{out}.txt',
            set=['transcript','gene'], out=['raw_counts', 'fpkm','tpm']
        )