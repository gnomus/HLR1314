#!/bin/sh

for n in 1 2 3
do
    for t in 1 2 3 4 5 6 7 8 9 10 11 12
    do
        sbatch -p cluster -o messung-threads_$t-run_$n.out job_script partdiff-posix $t 2 512 2 2 1024
    done
    sbatch -p cluster -o messung-seq-run$n.out job_script partdiff-seq 1 2 512 2 2 1024
    sbatch -p magny -o messung-threads_48-run_$n.out job_script partdiff-posix 48 2 512 2 2 1024
done
