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

sn  flowcell    lane    file
s1  fc1 l1  path/to/l1.fq.gz
s1  fc1 l3  path/to/l3.fq.gz
...

With a separate line for each lane.

See the local example file "chm13.mapping.manifest" for details.
"""

import pandas as pd

__author__ = ""
__license__ = "MIT"

configfile: "config.json"

manifest = pd.read_csv(config["manifest"], header=0, sep="\t")

shell.prefix("source config.sh; ")

def lanes_from_sample(wildcards):
    manifest_sn = manifest.loc[manifest.sn == wildcards.sample]
    bams = "mapping/" + wildcards.reference + "/" + manifest_sn.sn + "/" + manifest_sn.flowcell + "/" + manifest_sn.lane + ".bam"
    return bams.unique().tolist()

def get_files(wildcards):
    return manifest.loc[(manifest.sn == wildcards.sample) & (manifest.flowcell == wildcards.flowcell) & (manifest.lane == wildcards.lane), "file"].tolist()

from snakemake.exceptions import MissingInputException

rule all:
    input: expand("mapping/{reference}/metrics/{sample}.{type}.txt", reference = config["references"], sample = manifest["sn"].tolist(), type = ["insert_size_metrics", "flagstat", "idxstats"]),
           expand("mapping/{reference}/merged/{sample}.bam.bai", reference = config["references"], sample = manifest["sn"].tolist())
    params: sge_opts = "-N do_isize"

rule get_flagstat:
    input: "mapping/{reference}/merged/{sample}.bam", "mapping/{reference}/merged/{sample}.bam.bai"
    output: "mapping/{reference}/metrics/{sample}.flagstat.txt"
    params: sge_opts = "-l mfree=8G -N flagstat"
    shell:
        "samtools flagstat {input[0]} > {output}"
   

rule get_idxstats:
    input: "mapping/{reference}/merged/{sample}.bam", "mapping/{reference}/merged/{sample}.bam.bai"
    output: "mapping/{reference}/metrics/{sample}.idxstats.txt"
    params: sge_opts = "-l mfree=8G -N idxstats"
    shell:
        "samtools idxstats {input[0]} > {output}"

rule collect_isize_metrics:
    input: "mapping/{reference}/merged/{sample}.bam"
    output: "mapping/{reference}/metrics/{sample}.insert_size_metrics.txt"
    params: sge_opts = "-N collect_isize -l mfree=8G"
    shell:
        """java -Xmx8G -jar $PICARD_DIR/CollectInsertSizeMetrics.jar I={input} O={output} H={output}.hist.pdf"""

rule index_merged_bams:
    input: "mapping/{reference}/merged/{sample}.bam"
    output: "mapping/{reference}/merged/{sample}.bam.bai"
    params: sge_opts="-l mfree=8G -N index_bam"
    shell: "samtools index {input}"

rule merge_bams:
    input: lanes_from_sample
    output: "mapping/{reference}/merged/{sample}.bam"
    params: sge_opts="-l mfree=4G -pe serial 8 -N merge_bam"
    shell:
        "samtools merge -@ 8 {output} {input}"

rule bwa_mem_map_and_mark_dups:
    input:  lambda wildcards: config["references"][wildcards.reference],
            get_files
    output:
        "mapping/{reference}/{sample}/{flowcell}/{lane}.bam"
    params:
        sample="{sample}",
        flowcell="{flowcell}",
        custom=config.get("params_bwa_mem", ""),
        sge_opts="-l mfree=4G -pe serial 12 -N bwa_mem_map",
        bwa_threads = "8",
        samtools_threads = "4", samtools_memory = "8G"
    log:
        "mapping/log/{reference}/{sample}/{flowcell}/{lane}.log"
    shell:
        "bwa mem {params.custom} "
        r"-R '@RG\tID:{params.flowcell}_{wildcards.lane}\t"
        r"SM:{params.sample}\tLB:{params.sample}\tPL:{config[platform]}\tPU:{params.flowcell}' "
        "-t {params.bwa_threads} {input} 2> {log} "
        "| samblaster | samtools sort -@ {params.samtools_threads} -m {params.samtools_memory} -O bam -T $TMPDIR/{wildcards.lane} -o {output}; samtools index {output}"
