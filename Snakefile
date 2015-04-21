# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

"""
Read mapping with BWA MEM (http://bio-bwa.sourceforge.net).

Expects a config.json file like this:

{
    "references": {
        "hg19": "/net/eichler/vol2/eee_shared/assemblies/hg19/indexes/bwa_0.7.5a-r405/ucsc.hg19.fasta"
    },
    "manifest": "/path/to/manifest/directory",
    "platform": "ILLUMINA",
    "params_bwa_mem": "-M"
}

params_bwa_mem is optional

Expects a manifest file in tab-delimited format like this:

sample  flowcell    lane    file
s1  fc1 l1  path/to/l1.fq.gz
s1  fc1 l3  path/to/l3.fq.gz
...

With a separate line for each lane.
"""

import pandas as pd

__author__ = ""
__license__ = "MIT"

configfile: "config.json"

manifest = pd.read_csv(config["manifest"], header=0, sep="\t")

def lanes_from_sample(wildcards):
    return map(lambda x: "mapping/%s/lanes/%s.bam" % (wildcards.reference, x), manifest.loc[manifest.sample == wildcards.sample, "lane"].tolist())

def get_files(wildcards):
    return manifest.loc[manifest.lane == wildcards.lane, "file"].tolist()

from snakemake.exceptions import MissingInputException

rule all:
    input: expand("mapping/{reference}/metrics/{sample}.insert_size_metrics.txt", reference = config["references"], sample = manifest["sample"].tolist()),
           expand("mapping/{reference}/merged/{sample}.bam.bai", reference = config["references"], sample = manifest["sample"].tolist())
    params: sge_opts = "-N do_isize"

rule collect_isize_metrics:
    input: "mapping/{reference}/merged/{sample}.bam"
    output: "mapping/{reference}/metrics/{sample}.insert_size_metrics.txt"
    params: sge_opts = "-N collect_isize -l mfree=8G"
    shell:
        """source config.sh; java -Xmx8G -jar $PICARD_DIR/CollectInsertSizeMetrics.jar I={input} O={output} H={output}.hist"""

rule index_merged_bams:
    input: "mapping/{reference}/merged/{sample}.bam"
    output: "mapping/{reference}/merged/{sample}.bam.bai"
    params: sge_opts="-l mfree=8G -N index_bam"
    shell: "source config.sh; samtools index {input}"

rule merge_bams:
    input: lanes_from_sample
    output: "mapping/{reference}/merged/{sample}.bam"
    params: sge_opts="-l mfree=4G -pe serial 8 -N merge_bam"
    shell:
        "source config.sh; samtools merge -@ 8 {output} {input}"

rule bwa_mem_map_and_mark_dups:
    input:  lambda wildcards: config["references"][wildcards.reference],
            get_files
    output:
        "mapping/{reference}/lanes/{lane}.bam"
    params:
        sample=lambda wildcards: manifest.loc[manifest.lane == wildcards.lane, "sample"].tolist()[0],
        flowcell=lambda wildcards: manifest.loc[manifest.lane == wildcards.lane, "flowcell"].tolist()[0],
        custom=config.get("params_bwa_mem", ""),
        sge_opts="-l mfree=4G -pe serial 12 -N bwa_mem_map",
        bwa_threads = "8",
        samtools_threads = "4", samtools_memory = "8G"
    log:
        "mapping/log/{reference}/lanes/{lane}.log"
    shell:
        "source config.sh; "
        "bwa mem {params.custom} "
        r"-R '@RG\tID:{params.flowcell}_{wildcards.lane}\t"
        r"SM:{params.sample}\tLB:{params.sample}\tPL:{config[platform]}\tPU:{params.flowcell}' "
        "-t {params.bwa_threads} {input} 2> {log} "
        "| samblaster | samtools sort -@ {params.samtools_threads} -m {params.samtools_memory} -O bam -T /var/tmp/{wildcards.lane} -o {output}; samtools index {output}"
