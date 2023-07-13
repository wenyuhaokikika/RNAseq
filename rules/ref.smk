'''
@File    :   ref.smk
@Time    :   2023/07/10 11:15:14
@Author  :   wenyuhao 
@Version :   1.0
@Contact :   925201392@qq.com
@Desc    :   None
'''

# here put the import lib
import os
from functools import partial
from snakemake.shell import shell

def refPath(*args):
    return os.path.join(config['reference']['refDir'],*args)

def workPath(*args):
    return os.path.join(config['workspace'],*args)

#https://www.biostars.org/p/339405/
# https://bioinformatics.stackexchange.com/questions/540/what-ensembl-genome-version-should-i-use-for-alignments-e-g-toplevel-fa-vs-p
# 一个soft mask和hard mask和无mask。soft mask是小写，hard mask是N，无mask是大写
rule genome:
    output:
        fasta = refPath("genome.fa.gz"),
        gtf = refPath("genome.gtf.gz")
    params:
        release = config['reference']['release'],
        version = config['reference']['version']
    cache: True
    priority: 10
    retries: 3
    log:
        workPath('log/ref/genome.log'),
    shell:
        '''
        wget https://ftp.ensembl.org/pub/release-{params.release}/gtf/homo_sapiens/Homo_sapiens.{params.version}.{params.release}.chr_patch_hapl_scaff.gtf.gz -O {output.gtf}  2> {log}
        wget https://ftp.ensembl.org/pub/release-{params.release}/fasta/homo_sapiens/dna/Homo_sapiens.{params.version}.dna_sm.toplevel.fa.gz -O {output.fasta}  2> {log}
        '''

rule sequenceDownload:
    output:
        cdna =  refPath('cdna.fa.gz'),
        pep = refPath('pep.fa.gz'),
        cds = refPath('cds.fa.gz')
    params:
        release = config['reference']['release'],
        version = config['reference']['version']
    cache: True
    priority: 10
    log:
        workPath('log/ref/sequenceDownload.log')
    retries: 3
    shell:
        '''
        wget https://ftp.ensembl.org/pub/release-{params.release}/fasta/homo_sapiens/cdna/Homo_sapiens.{params.version}.cdna.all.fa.gz -O {output.cdna}   2> {log}
        wget https://ftp.ensembl.org/pub/release-{params.release}/fasta/homo_sapiens/pep/Homo_sapiens.{params.version}.pep.all.fa.gz -O {output.pep} 2> {log}
        wget https://ftp.ensembl.org/pub/release-{params.release}/fasta/homo_sapiens/cds/Homo_sapiens.{params.version}.cds.all.fa.gz -O  {output.cds} 2> {log}
        '''
checkpoint STARindex:
    output:
        directory(refPath('STARindex'))
    input:
        fasta = rules.genome.output.fasta, #refPath("genome.fa.gz")
        gtf = rules.genome.output.gtf #refPath("genome.gtf.gz")
    cache: True
    threads:
        threads = config['threads']
    priority: 50
    log:
        workPath('log/ref/STARindex.log')
    retries: 3
    shell:
        'STAR'
        '--runMode genomeGenerate'
        '--genomeDir {output}'
        '--runThreadN {threads}'
        '--genomeFastaFiles {input.fasta}'
        '--sjdbGTFfile {input.gtf}'
        '--sjdbOverhang 125'
        ' 2> {log}'

rule annotation:
    output:
        directory(refPath('annotation'))
    params:
        release = config['reference']['release'],
        version = config['reference']['version']
    cache: True
    priority: 50
    log:
        workPath('log/ref/annotation.log')
    retries: 3
    shell:
        '''
        wget https://ftp.ensembl.org/pub/release-{params.release}/tsv/homo_sapiens/Homo_sapiens.GRCh38.{params.release}.uniprot.tsv.gz -P {output} 2> {log}
        wget https://ftp.ensembl.org/pub/release-{params.release}/tsv/homo_sapiens/Homo_sapiens.GRCh38.{params.release}.karyotype.tsv.gz -P  {output} 2> {log}
        wget https://ftp.ensembl.org/pub/release-{params.release}/tsv/homo_sapiens/Homo_sapiens.GRCh38.{params.release}.entrez.tsv.gz -P {output} 2> {log}
        wget https://ftp.ensembl.org/pub/release-{params.release}/tsv/homo_sapiens/Homo_sapiens.GRCh38.{params.release}.ena.tsv.gz -P {output} 2> {log}
        wget https://ftp.ensembl.org/pub/release-{params.release}/tsv/homo_sapiens/Homo_sapiens.GRCh38.{params.release}.canonical.tsv.gz -P {output} 2> {log}
        wget https://ftp.ensembl.org/pub/release-{params.release}/tsv/homo_sapiens/README_ENA.tsv -P  {output} 2> {log}
        wget https://ftp.ensembl.org/pub/release-{params.release}/tsv/homo_sapiens/README_refseq.tsv -P  {output} 2> {log}
        wget https://ftp.ensembl.org/pub/release-{params.release}/tsv/homo_sapiens/README_entrez.tsv -P  {output} 2> {log}
        wget https://ftp.ensembl.org/pub/release-{params.release}/tsv/homo_sapiens/README_uniprot.tsv -P  {output} 2> {log}
        '''