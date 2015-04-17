# bwa_mem_mapping
Fastq to bam mapping pipeline with marked duplicates, indexing, and insert size metric collection

This is a Snakemake-based pipeline for mapping with BWA mem. The pipeline automatically marks duplicates, indexes bams, and calculates insert size metrics, all of which are useful for downstream analysis.
The Snakefile is derived from johanneskoester's example here:
https://bitbucket.org/johanneskoester/snakemake-workflows/src/e03c8a9fb7f256ed498b0a5984e9ae738a8f22e1/bio/ngs/rules/mapping/bwa_mem.rules?at=master
