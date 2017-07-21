#!/usr/bin/sh
#shell script for submitting VAMF jobs to BSUB system

max=3
for (( i=1; i <= $max; ++i ))
do
    echo "data/MGH26/batch_separability/subsets/${i}.tsv"
done

#which Rscript
#Rscript --version
fname="test"
prefix="data/MGH26/batch_separability"
mkdir -p ${prefix}/bsub_out
bsub -q short -o ${prefix}/bsub_out/${fname}.out -e ${prefix}/bsub_out/${fname}.err -J VAMF${fname} Rscript ../algs/vamf_script.R -i ${prefix}/subsets/${fname}.tsv -o ${prefix}/vamf_out/${fname}.tsv
