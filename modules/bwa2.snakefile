
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
    threads: 16
    params:
        tmp="analysis/bwa/{sample}/temp",
        ref_genome=config['ref_genome']
    shell:
        "mkdir -p {params.tmp} && "
        "sambamba sort" 
        " --tmpdir {params.tmp}"
        " --out {output.bam}"
        " --show-progress"
        " --nthreads {threads}"
        " {input.bam}"

rule dedupbam:
    input:
        bam="analysis/bwa/{sample}/{sample}_{lane}.sort.bam"
    output:
        bam="analysis/bwa/{sample}/{sample}_{lane}.dedup.bam",
    wildcard_constraints:
        sample="[^/]+",
        lane="[^/]+"
    benchmark: "benchmarks/dedupbam.{sample}_{lane}.txt"
    threads: 16
    params:
        tmp="analysis/bwa/{sample}/temp",
        java_opts="-Xms50G"
    shell:
        "sambamba markdup"
        " --remove-duplicates"
        " --nthreads {threads}"
        " --show-progress"
        " --tmpdir {params.tmp}"
        " {input.bam}"
        " {output.bam}"


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
    return " ".join(inputs)

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
        inputs=mergebam_input_args,
    shell:
        "sambamba merge"
        " --nthreads {threads}"
        " --show-progress"
        " {output.bam}"
        " {params.inputs} && "
        "sambamba index"
        " --nthreads {threads}"
        " --show-progress"
        " {output.bam}"
        " {output.bai}"


