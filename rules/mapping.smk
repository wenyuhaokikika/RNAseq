'''
@File    :   mapping.smk
@Time    :   2023/01/30 14:21:50
@Author  :   wenyuhao 
@Version :   1.0
@Contact :   925201392@qq.com
@Desc    :   None
'''

# here put the import lib

import os
rule mapping:
    output:
        config['workspace'] + '/samples/{sample}/mapping/{library}_Aligned.out.bam'
    input:
        fq1=rules.trim.output.fq1,
        fq2=rules.trim.output.fq2
    log:
        config['workspace'] + '/log/mapping/{sample}/{library}_star.log'
    params:
        genome_dir=config['genome']['star_index'],
        prefix=config['workspace'] + '/samples/{sample}/mapping/{library}_'
    threads:
        20 if workflow.cores > 20 else workflow.cores
    priority:
        10
    shell:
        'STAR --readFilesIn {input.fq1} {input.fq2}'
        ' --outSAMattrRGline ID:{wildcards.library} SM:{wildcards.sample} LB:{wildcards.library} PL:ILLUMINA'
        ' --alignIntronMax 1000000'
        ' --alignIntronMin 20'
        ' --alignMatesGapMax 1000000'
        ' --alignSJDBoverhangMin 1'
        ' --alignSJoverhangMin 8'
        ' --alignSoftClipAtReferenceEnds Yes'
        ' --chimJunctionOverhangMin 15'
        ' --chimMainSegmentMultNmax 1'
        ' --chimOutType Junctions SeparateSAMold WithinBAM SoftClip'
        ' --chimSegmentMin 15'
        ' --genomeDir {params.genome_dir}'
        ' --genomeLoad NoSharedMemory'
        ' --limitSjdbInsertNsj 1200000'
        ' --outFileNamePrefix {params.prefix}'
        ' --outFilterIntronMotifs None'
        ' --outFilterMatchNminOverLread 0.33'
        ' --outFilterMismatchNmax 999'
        ' --outFilterMismatchNoverLmax 0.1'
        ' --outFilterMultimapNmax 20'
        ' --outFilterScoreMinOverLread 0.33'
        ' --outFilterType BySJout'
        ' --outSAMattributes NH HI AS nM NM ch'
        ' --outSAMstrandField intronMotif'
        ' --outSAMtype BAM Unsorted'
        ' --outSAMunmapped Within'
        ' --quantMode TranscriptomeSAM GeneCounts'
        ' --readFilesCommand zcat'
        ' --runThreadN {threads}'
        ' --twopassMode Basic >{log} 2>&1'
