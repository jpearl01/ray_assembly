#!/bin/bash
prefix=$1
java -jar /opt/Trimmomatic-0.33/trimmomatic-0.33.jar \
PE \
-threads 16 \
-phred33 \
rawdata/$prefix.1.fq.gz \
rawdata/$prefix.2.fq.gz \
cleandata/$prefix.1.fq \
cleandata/$prefix.unpaired.1.fq \
cleandata/$prefix.2.fq \
cleandata/$prefix.unpaired.2.fq \
ILLUMINACLIP:IlluminaContaminants.fa:2:30:10 \
1>log/trim.$prefix.stdout \
2>log/trim.$prefix.stderr \

