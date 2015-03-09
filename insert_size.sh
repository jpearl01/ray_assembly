#!/bin/bash
prefix=$1
java -Xmx10g -jar bin/picard-tools-1.97/CollectInsertSizeMetrics.jar \
VALIDATION_STRINGENCY=SILENT \
HISTOGRAM_FILE="insert_size/insert_size.$prefix.pdf" \
INPUT="bam/$prefix.sort.bam" \
OUTPUT="insert_size/insert_size.$prefix.txt" \
ASSUME_SORTED=false \
1>log/insert_size.metrics.$prefix.stdout \
2>log/insert_size.metrics.$prefix.stderr

