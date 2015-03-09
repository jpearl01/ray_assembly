#!/bin/bash
prefix=$1
bin/cope-v1.1.2/cope/cope \
-a cleandata/$prefix.1.fq \
-b cleandata/$prefix.2.fq \
-o cleandata/$prefix.connected.fq \
-2 cleandata/$prefix.unconnected.1.fq \
-3 cleandata/$prefix.unconnected.2.fq \
-l 10 \
-u 100 \
-s 33 \
-m 0 \
1>log/connect.$prefix.stdout \
2>log/connect.$prefix.stderr


#Program:	COPE (Connecting Overlapped Pair-End reads)
#AUTHOR:		BGI shenzhen
#Version:	COPE v1.1.0
#Compile Date:	Sep 17 2012 time: 14:26:36
#Contact: 	liubinghang@genomics.org.cn; yuanjianying@genomics.org.cn
#Usage:	cope [option]
#	-a  <str>   query read1 file, (*.fq, *.fa, *.fq.gz, *.fa.gz)
#	-b  <str>   query read2 file
#	-o  <str>   output connected file in *.fq or *fa
#	-2  <str>   output fail connected read1.fq
#	-3  <str>   output fail connected read2.fq
#	-l  <int>   lower bound overlap length, 8 or 10 is preferred for 100bp read with 170~180 insert size, default=10
#	-u  <int>   higher bound overlap length, 70 is preferred for 100 bp reads, default=70
#	-c  <float>   match ratio cutoff, default=0.75
#	-d  <float>   match max2_ratio/max_ratio, important for mode 0, default=0.7
#	-B  <float>   ratio cut off value of Base-quality=2, default=0.9
#	-N  <int>   N filter( 1 filte 0 not filter), default=1
#	-T  <int>   read pair number threshold for parameter training or testing: 0
#	-s  <int>   the smallest ASCII value used to represent base quality in fastq files. 33 for the later Illumina platforms and Sanger reads, and 64 for the earlier Illumina platforms. default=64
#	-m  <int>   mode type: 0 simple connect; 1 k-mer frequency assisted connection; 2 Auxiliary reads and cross connection for left reads from mode 1; 3 full mode ,default=3
#
#	when set mode large than 0, set the following options:
#	-k  <int>   Kmer size, default=17
#	-t  <str>   compressed kmer frequency table
#	-f  <str>   kmer frequency table len, for parallelly compressed kmer freq table
#	-j  <str>   compressed kmer frequency table counted by jellyfish
#	-L  <int>   set the threshold of freq-low(Kmer with frequency lower than(<=) this threshold will not be considered for base selection), default=3
#	-M  <int>   set the threshold of freq-normal(Kmer whith frequency lower than(<=) this threshold will not be considered for spanning Kmer selection), default=10
#	-H  <int>   set the threshold of freq-high(Kmer whith frequency higher than(>=) this threshold will not be considered for spanning Kmer selection, 2 times the average kmer frequency is preferred), default=60
#	-R  <int>   set the threshold of spanning k-mers number of extremely high frequency, threshold=3 is preferred for high repeat genome, set -1 for disable this threshold, default=-1
#	-K  <str>   input pair-kmer file for mode=2, default=./pair.kmer.list.gz
#
#	when set mode=2 or 3, you also need to set the following options:
#	-D  <int>   set the batch number of cross-information file(cross.read.inf*.gz) output, default=1
#	-r  <int>   use connected reads(1) or raw reads(0) to find cross reads, default=1
#	-F  <str>   input extra read-files list(each line for one read file) while need to find cross reads with the help of other read-files
#
#	-h          print help information
#
#Example:
# ##simple connect mode: ./cope -a read1.fq -b read2.fq -o connect.fa -2 left1.fq -3 left2.fq -m 0 >cope.log 2>cope.error
# ##k-mer frequency assisted connection: ./cope -a read1.fq -b read2.fq -o connect.fq -2 left1.fq -3 left2.fq -m 1 -t kmer_table.cz -f kmer_table.cz.len >cope.log 2>cope.error
# ##Auxiliary reads and cross connection for left reads from mode 1: ./cope -a read1.fq -b read2.fq -o connect.fq -2 left1.fq -3 left2.fq -m 2 -t kmer_table.cz -f kmer_table.cz.len >cope.log 2>cope.error
# ##connect by full mode: ./cope -a read1.fq -b read2.fq -o connect.fq -2 left1.fq -3 left2.fq -m 3 -t kmer_table.cz -f kmer_table.cz.len >cope.log 2>cope.error
#
#
