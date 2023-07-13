'''
@File    :   common.smk
@Time    :   2023/07/06 15:17:24
@Author  :   wenyuhao 
@Version :   1.0
@Contact :   925201392@qq.com
@Desc    :   None
'''

# here put the import lib
import os
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

configfile: "config/config.yaml"

# contigs in reference genome

# units = pd.read_table(config["units"], dtype=str,comment='#').set_index(
#     ["sample", "unit"], drop=False
# )

# units.index = units.index.set_levels(
#     [i.astype(str) for i in units.index.levels]
# )  # enforce str in index
units = pd.read_table(config["units"], dtype=str,comment='#')
samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
libraryD = units.set_index('sample')['library'].to_dict()

wildcard_constraints:
    sample="|".join(samples.index),
    unit="|".join(units["unit"]),

def get_contigs():
    with checkpoints.genome_faidx.get().output[0].open() as fai:
        return pd.read_table(fai, header=None, usecols=[0], squeeze=True, dtype=str)

def get_fastq(wildcards):
    """
    Get fastq files of given sample-unit.
    need columns ["sample",'unit',"fq1","fq2"].
    """
    # fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    fastqs = units.loc[(units["sample"]==wildcards.sample)&(units["unit"]==wildcards.unit),:].iloc[0]#["fq1", "fq2"].dropna()
    #print(fastqs,dict(fastqs))
    if not pd.isna(fastqs.fq2):#len(fastqs) == 2:
        return {"fq1": fastqs.fq1, "fq2": fastqs.fq2}
    return {"fq1": fastqs.fq1}

def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(
            "results/trimmed/{sample}-{unit}.{group}.fastq.gz",
            group=[1, 2],
            **wildcards
        )
    # single end sample
    return "results/trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        "results/recal/{sample}-{unit}.bam",
        sample=wildcards.sample,
        unit=units.loc[wildcards.sample].unit,
    )
