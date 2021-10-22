
wkdir="/data08/project/dev_test_fullgenome_v2"
conffn="/data08/project/dev_test_fullgenome_v2/dev.config.yaml"
queue="kdna.q"

mkdir -p ${wkdir}/logs
cmd="snakemake \
    --cores all \
    --printshellcmds \
    --snakefile /TBI/Share/BioPeople/siyoo/Pipelines/kdna_qc_pipeline/Snakemake \
    --config 'ref=/TBI/Share/BioPeople/siyoo/Pipelines/kdna_qc_pipeline/refs/ref_kdna.yaml' \
    --cluster-config /TBI/Share/BioPeople/siyoo/Pipelines/kdna_qc_pipeline/envs/cluster.json \
    --max-jobs-per-second 5 \
    --jobs 5 \
    --configfile ${conffn} \
    --directory ${wkdir} \
    --cluster 'qsub -V -o ${wkdir}/{cluster.output} -e ${wkdir}/{cluster.error} -pe smp {cluster.threads} -N {cluster.jobName} -S /bin/bash -q ${queue}'"
echo $cmd


