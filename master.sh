#!/bin/bash
#bin/bwa-0.7.5a/bwa index ref/NP.fa

prefix=$1
./trim.sh $prefix
./align_to_ref.sh $prefix
./insert_size.sh $prefix
./connect.sh $prefix
./error_correct.sh $prefix
./kmergenie.sh $prefix
./ray_assemble.sh $prefix
