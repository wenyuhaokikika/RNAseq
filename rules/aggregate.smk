'''
@File    :   aggregate.smk
@Time    :   2023/01/30 14:16:31
@Author  :   wenyuhao 
@Version :   1.0
@Contact :   925201392@qq.com
@Desc    :   None
'''

def find_all_library(wildcards):
    return [
        config['workspace'] + f'/samples/{wildcards.sample}/mapping/{unit}_Aligned.out.bam'
        for unit in units['unit']
    ]

rule countGene:
    output:
        config['workspace'] + '/aggregate/{sample}_gene_raw_count.txt'
    input:
        find_all_library
    log:
        config['workspace'] + '/log/aggregate/{sample}/{sample}_gene_htseq_count.log'
    params:
        gtf=rules.genome.output.gtf #config['genome']['gtf']
    shell:
        'htseq-count -f bam -r name -s no -a 10 -t exon -i gene_id -m intersection-nonempty'
        ' {input} {params.gtf} > {output} 2>{log}'

rule countTranscript:
    output:
        config['workspace'] + '/aggregate/{sample}_transcipt_raw_count.txt'
    input:
        find_all_library
    log:
        config['workspace'] + '/log/aggregate/{sample}/{sample}_transcript_htseq_count.log'
    params:
        gtf=rules.genome.output.gtf #gtf=config['genome']['gtf']
    shell:
        'htseq-count -f bam -r name -s no -a 10 -t exon -i transcript_id -m intersection-nonempty'
        ' {input} {params.gtf} > {output} 2>{log}'

rule merge_count:
    output:
        rawGene=config['workspace'] + '/aggregate/all_sample_gene_raw_counts.txt',
        fpkmGene=config['workspace'] + '/aggregate/all_sample_gene_fpkm.txt',
        tpmGene=config['workspace'] + '/aggregate/all_sample_gene_tpm.txt',
        rawTranscript=config['workspace'] + '/aggregate/all_sample_transcript_raw_counts.txt',
        fpkmTranscript=config['workspace'] + '/aggregate/all_sample_transcript_fpkm.txt',
        tpmTranscript=config['workspace'] + '/aggregate/all_sample_transcript_tpm.txt',
    input:
        expand(
            config['workspace'] + '/aggregate/{sample}_{biounit}_raw_count.txt',
            sample=samples['sample'], #config['samples']
            biounit=['gene','transcipt']
        )
    params:
        names=expand('{sample}', sample=samples['sample']),
        gtf = rules.genome.output.gtf
    run:
        # https://github.com/reneshbedre/bioinfokit/blob/master/bioinfokit/analys.py
        import pandas as pd
        import numpy as np
        import gtfparse
        import qnorm
        from typing import Tuple
        def parse_name(name)->Tuple[str,str,str]:
            r'''return (sample,biounit,path)'''
            return (*name.split('/')[-1].split('_')[:2],name)
        paths = [parse_name(i) for i in input]
        geneCounts = pd.concat([pd.read_table(i[-1], names=['id', i[0]], index_col=0) for i in paths if i[1] == 'gene'], axis=1).drop('id',axis=1)
        transCounts = pd.concat([pd.read_table(i[-1], names=['id', i[0]], index_col=0) for i in paths if i[1] == 'transcipt'], axis=1).drop('id',axis=1)
        geneCounts = geneCounts.loc[~geneCounts.index.str.startswith('__')]
        transCounts = transCounts.loc[~transCounts.index.str.startswith('__')]
        geneCounts.to_csv(output.rawGene, sep='\t', index_label='')
        transCounts.to_csv(output.rawTranscript, sep='\t', index_label='')

        gtf = gtfparse.read_gtf(params.gtf) #config['genome']['gtf'])
        transcripts = gtf.loc[
            (gtf['feature'] == 'transcript'), 'transcript_id'
        ].drop_duplicates().values
        genes = gtf.loc[
            (gtf['feature'] == 'transcript'), 'gene_id'
        ].drop_duplicates().values

        totalTrans = transCounts.loc[transcripts].sum(0)
        totalGene = geneCounts.loc[genes].sum(0)

        def exon_length(exon):
            exon = exon.sort_values(['start', 'end'])
            exon['group'] = (exon['start'] > exon['end'].shift()).cumsum()
            width = exon.groupby('group').agg({'start': 'min', 'end': 'max'})
            return (width['end'] - width['start']).sum()

        exons = gtf[gtf['feature'] == 'exon']
        lengthTrans = exons.groupby('transcript_id').apply(exon_length)
        lengthGene = exons.groupby('gene_id').apply(exon_length)
        (transCounts.loc[transcripts,:] * 10 ** 9 / totalTrans).div(lengthTrans[transcripts], axis=0).to_csv(output.fpkmTranscript, sep='\t', index_label='')
        (geneCounts.loc[genes,:] * 10 ** 9 / totalGene).div(lengthGene[genes], axis=0).to_csv(output.fpkmGene, sep='\t', index_label='')

        def tpm(counts,length,index):
            a = counts.loc[index,:].div(length[index], axis=0) * 1e3
            return (a * 1e6) / a.sum()

        tpm(transCounts, lengthTrans,transcripts).to_csv(output.tpmTranscript, sep='\t', index_label='')
        tpm(geneCounts, lengthGene,genes).to_csv(output.tpmGene, sep='\t', index_label='')
