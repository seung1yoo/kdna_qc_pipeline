
rule haplotypecaller:
    input:
        bam="analysis/bwa/{sample}/{sample}.merge.bam"
    output:
        vcf="analysis/haplotypecaller/{sample}/{sample}.vcf.gz"
    wildcard_constraints:
        sample="[^/]+",
    benchmark: "benchmarks/haplotypecaller.{sample}.txt"
    threads: 16
    params:
        tmp="analysis/haplotypecaller/{sample}/temp",
        java_opts="-Xmx4g"
    shell:
        "mkdir -p {params.tmp} && "
        "gatk"
        " --java-options '{params.java_opts}'"
        " HaplotypeCaller"
        " --tmp-dir {params.tmp}"
        " -R {config[ref_genome]}"
        " -I {input.bam}"
        " -O {output.vcf}"
        " --intervals {config[ref_region]}"
        " --dbsnp {config[ref_dbsnp]}"
        " --native-pair-hmm-threads {threads}"

