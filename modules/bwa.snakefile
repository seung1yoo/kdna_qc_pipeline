

rule bwa:
    input:
        fq1=lambda wildcards: config["samples"][wildcards.sample][wildcards.lane]["fq1"],
        fq2=lambda wildcards: config["samples"][wildcards.sample][wildcards.lane]["fq2"]
    output:
        bam="analysis/bwa/{sample}/{sample}_{lane}.bam"
    wildcard_constraints:
        sample="[^/]+",
        lane="[^/]+"
    benchmark: "benchmarks/bwa.{sample}_{lane}.txt"
    threads: 16
    params:
        tmp="analysis/bwa/{sample}/temp",
        seed="45",
        rg="\'@RG\\tID:{sample}.{lane}\\tPL:illumina\\tSM:{sample}.{lane}\'",
    shell:
        "mkdir -p {params.tmp} && "
        "bwa mem -M"
        " -R {params.rg}"
        " -t {threads}"
        " -k {params.seed}"
        " {config[ref_bwaidx]}"
        " {input.fq1} {input.fq2} | "
        "samtools view -Sb -F 0x100 > {output.bam}"

rule sortbam:
    input:
        bam="analysis/bwa/{sample}/{sample}_{lane}.bam"
    output:
        bam="analysis/bwa/{sample}/{sample}_{lane}.sort.bam"
    wildcard_constraints:
        sample="[^/]+",
        lane="[^/]+"
    benchmark: "benchmarks/sortbam.{sample}_{lane}.txt"
    threads: 1
    params:
        tmp="analysis/bwa/{sample}/temp",
        ref_genome=config['ref_genome']
    shell:
        "mkdir -p {params.tmp} && "
        "gatk SortSam"
        " --TMP_DIR {params.tmp}"
        " --INPUT {input.bam}"
        " --OUTPUT {output.bam}"
        " --SORT_ORDER coordinate"
        " --CREATE_INDEX true"
        " --CREATE_MD5_FILE true"


rule dedupbam:
    input:
        bam="analysis/bwa/{sample}/{sample}_{lane}.sort.bam"
    output:
        bam="analysis/bwa/{sample}/{sample}_{lane}.dedup.bam",
        bai="analysis/bwa/{sample}/{sample}_{lane}.dedup.bai",
        txt="analysis/bwa/{sample}/{sample}_{lane}.dedup.txt"
    wildcard_constraints:
        sample="[^/]+",
        lane="[^/]+"
    benchmark: "benchmarks/dedupbam.{sample}_{lane}.txt"
    threads: 1
    params:
        tmp="analysis/bwa/{sample}/temp",
        java_opts="-Xms50G"
    shell:
        "mkdir -p {params.tmp} && "
        "gatk"
        " --java-options '{params.java_opts}'"
        " MarkDuplicates"
        " --TMP_DIR {params.tmp}"
        " --INPUT {input.bam}"
        " --OUTPUT {output.bam}"
        " --METRICS_FILE {output.txt}"
        " --VALIDATION_STRINGENCY SILENT"
        " --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500"
        " --REMOVE_SEQUENCING_DUPLICATES true"
        " --CREATE_INDEX true"
        " --CREATE_MD5_FILE true"


def mergebam_input(wildcards):
    inputs = list()
    sample = wildcards.sample
    for lane in config["samples"][wildcards.sample]:
        fn = f"analysis/bwa/{sample}/{sample}_{lane}.dedup.bam"
        inputs.append(fn)
    return inputs

def mergebam_input_args(wildcards):
    inputs = list()
    sample = wildcards.sample
    for lane in config["samples"][wildcards.sample]:
        fn = f"analysis/bwa/{sample}/{sample}_{lane}.dedup.bam"
        inputs.append(fn)
    print(inputs)
    return " --INPUT ".join(inputs)


rule mergebam:
    input:
        bam=mergebam_input
    output:
        bam="analysis/bwa/{sample}/{sample}.merge.bam",
        bai="analysis/bwa/{sample}/{sample}.merge.bai"
    wildcard_constraints:
        sample="[^/]+",
    benchmark: "benchmarks/mergebam.{sample}.txt"
    threads: 16
    params:
        tmp="analysis/bwa/{sample}/temp",
        java_opts="-Xms50G",
        inputs=mergebam_input_args
    shell:
        "mkdir -p {params.tmp} && "
        "gatk"
        " --java-options '{params.java_opts}'"
        " MergeSamFiles"
        " --TMP_DIR {params.tmp}"
        " --INPUT {params.inputs}"
        " --OUTPUT {output.bam}"
        " --SORT_ORDER coordinate"
        " --USE_THREADING true"
        " --CREATE_INDEX true"

