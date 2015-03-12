#!/usr/bin/env bash

#Align the reads to the reference
/opt/bwa-0.7.10/bwa mem -t 16 -M ref/NP.fa $clean/$prefix.1.fq $clean/$prefix.2.fq 2>$log/align_to_ref.$prefix.stderr \
| bin/samtools-0.1.19/samtools view -Su - 2> $log/samtools.view.$prefix.stderr \
| bin/samtools-0.1.19/samtools sort -m 10G -l 0 -@ 8 - $bam/$prefix.sort 2>$log/samtools.sort.$prefix.stderr

#Insert size computation
java -Xmx10g -jar /opt/picard-tools-1.119/CollectInsertSizeMetrics.jar \
VALIDATION_STRINGENCY=SILENT \
HISTOGRAM_FILE="$is/insert_size.$prefix.pdf" \
INPUT="$bam/$prefix.sort.bam" \
OUTPUT="$is/insert_size.$prefix.txt" \
ASSUME_SORTED=false \
1>$log/insert_size.metrics.$prefix.stdout \
2>$log/insert_size.metrics.$prefix.stderr


#Error Correct the reads
mean=`grep -A 1 MEDIAN_INSERT_SIZE $is/insert_size.$prefix.txt | tail -n 1 | awk '{print $5-202}' | cut -f 1 -d \.`
stddev=`grep -A 1 MEDIAN_INSERT_SIZE $is/insert_size.$prefix.txt | tail -n 1 | awk '{print $6}' | cut -f 1 -d \.`

