'''
@File    :   Snakefile
@Time    :   2023/01/30 14:23:11
@Author  :   wenyuhao 
@Version :   1.0
@Contact :   925201392@qq.com
@Desc    :   None
'''

# here put the import lib
include: 'rules/preprocess.smk'
include: 'rules/mapping.smk'
include: 'rules/aggregate.smk'
include: 'rules/comparison.smk'


import os

BASE_DIR = os.path.dirname(workflow.snakefile)


rule all:
    input:
        expand(
            config['workspace'] + '/aggregate/all_sample_{set}_{out}.txt',
            set=['all', 'ccds'], out=['raw_counts', 'fpkm','tpm', 'tpm_qnorm']
        ),
        expand(
            config['workspace'] + '/aggregate/all_sample_{set}_pcaplot.{fmt}',
            fmt=config['plot_formats'], set=['all', 'ccds']
        ),
        expand(
            config['workspace'] + '/comparisons/{comparison}/{comparison}_{set}_result.txt',
            comparison=config['comparisons'], set=['all', 'ccds']
        ),
        expand(
            config['workspace'] + '/comparisons/{comparison}/{comparison}_gsea_{geneset}',
            comparison=config['comparisons'], geneset=config['genome']['geneset']
        )
