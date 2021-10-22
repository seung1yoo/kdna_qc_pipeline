


rule bam_stats:
    input:
        bam="analysis/bwa/{sample}/{sample}.merge.bam",
    output:
        stats="analysis/result/{sample}/{sample}.align.stats",
    wildcard_constraints:
        sample="[^/]+",
    benchmark: "benchmarks/bam_stats.{sample}.txt"
    threads: 16
    shell:
        "samtools stats"
        " --threads {threads}"
        " {input.bam} | "
        "grep ^SN | cut -f 2-"
        " > {output.stats}"



