#!/bin/bash
#bin/bwa-0.7.5a/bwa index ref/NP.fa

echo "This script takes two paired-end .tar.gz files and the mean and standard deviation of insert size"
echo "and creates a de novo assembly of the reads.  Note the files should begin with the same prefix"
echo "Usage:"
echo "./master.sh paired_file1 paired_file2 mean_insert stdev_insert"
echo "Example: ./master.sh paired.1.fq.gz paired.2.fq.gz 200 10"

filename1=$1
filename2=$2

full_path_fn1=$(readlink -f filename1)
full_path_fn2=$(readlink -f filename2)

prefix=$(basename $filename1)
extension="${filename1##*.}"
prefix="${prefix%.*}"

mean=$3
stdev=$4

clean=$prefix/cleandata
raw=$prefix/raw
log=$prefix/log
bam=$prefix/bam
is=$prefix/insert_size
ec=$prefix/ec_data
kg=$prefix/kmergenie
asmbly=$prefix/assembly

[ -d $raw ]||mkdir -p $raw
[ -d $clean ]||mkdir $clean
[ -d $log ]||mkdir $log
[ -d $bam ]||mkdir $bam
[ -d $is ]||mkdir $is
[ -d $ec ]||mkdir $ec
[ -d $kg ]||mkdir $kg
[ -d $asmbly ]||mkdir $asmbly

#Link the files into the new rawdata directory:
ln -s $full_path_fn1 $raw/$prefix.1.fq.gz
ln -s $full_path_fn2 $raw/$prefix.2.fq.gz


#Trim the reads in the input file
java -jar /opt/Trimmomatic-0.33/trimmomatic-0.33.jar \
PE \
-threads 16 \
-phred33 \
$raw/$prefix.1.fq.gz \
$raw/$prefix.2.fq.gz \
$clean/$prefix.1.fq \
$clean/$prefix.unpaired.1.fq \
$clean/$prefix.2.fq \
$clean/$prefix.unpaired.2.fq \
ILLUMINACLIP:IlluminaContaminants.fa:2:30:10 \
1>$log/trim.$prefix.stdout \
2>$log/trim.$prefix.stderr \

:'
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
'

#Connect paired-end reads
/opt/cope-v1.1.2/cope/cope \
-a $clean/$prefix.1.fq \
-b $clean/$prefix.2.fq \
-o $clean/$prefix.connected.fq \
-2 $clean/$prefix.unconnected.1.fq \
-3 $clean/$prefix.unconnected.2.fq \
-l 10 \
-u 100 \
-s 33 \
-m 0 \
1>$log/connect.$prefix.stdout \
2>$log/connect.$prefix.stderr

:'
#Error Correct the reads
mean=`grep -A 1 MEDIAN_INSERT_SIZE $is/insert_size.$prefix.txt | tail -n 1 | awk '{print $5-202}' | cut -f 1 -d \.`
stddev=`grep -A 1 MEDIAN_INSERT_SIZE $is/insert_size.$prefix.txt | tail -n 1 | awk '{print $6}' | cut -f 1 -d \.`
'

/opt/allpathslg/bin/ErrorCorrectReads.pl \
PHRED_ENCODING=33 \
PAIRED_READS_A_IN="$clean/$prefix.unconnected.1.fq" \
PAIRED_READS_B_IN="$clean/$prefix.unconnected.2.fq" \
PAIRED_SEP="$mean" \
PAIRED_STDEV="$stddev" \
PLOIDY=1 \
MAX_MEMORY_GB=200 \
THREADS=32 \
UNPAIRED_READS_IN="$clean/$prefix.connected.fq" \
READS_OUT="$ec/$prefix" \
1>$log/error_correct.$prefix.stdout \
2>$log/error_correct.$prefix.stderr



#Compute kmer info
ls $ec/$prefix.*.fastq >$kg/$prefix.reads.list

/opt/kmergenie-1.6950/kmergenie \
$kg/$prefix.reads.list \
-t 32 \
-o $kg/$prefix \
1>$log/kmergenie.$prefix.stdout \
2>$log/kmergenie.$prefix.stderr


#Assemble Data with ray
bestk=`tail -n 1 $log/kmergenie.$prefix.stdout | sed -e "s/best k: //"`

mpiexec -n 16 /opt/Ray-2.3.1/Ray \
-k $bestk \
-p $ec/$prefix.paired.A.fastq $ec/$prefix.paired.B.fastq \
-s $ec/$prefix.unpaired.fastq \
-o $asmbly/ray.k.$bestk.$prefix




