#!/bin/bash
prefix=$1
bin/bwa-0.7.5a/bwa mem -t 32 -M ref/NP.fa cleandata/$prefix.1.fq cleandata/$prefix.2.fq 2>log/align_to_ref.$prefix.stderr \
| bin/samtools-0.1.19/samtools view -Su - 2> log/samtools.view.$prefix.stderr \
| bin/samtools-0.1.19/samtools sort -m 10G -l 0 -@ 8 - bam/$prefix.sort 2>log/samtools.sort.$prefix.stderr

