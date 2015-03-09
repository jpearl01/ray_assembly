#!/bin/bash
prefix=$1

ls ec_data/$prefix.*.fastq >kmergenie/$prefix.reads.list

bin/kmergenie-1.5658/kmergenie \
kmergenie/$prefix.reads.list \
-t 32 \
-o kmergenie/$prefix \
1>log/kmergenie.$prefix.stdout \
2>log/kmergenie.$prefix.stderr

#KmerGenie
#
#Usage:
#    kmergenie <read_file> [options]
#
#Options:
#    --diploid    use the diploid model
#    --one-pass   skip the second pass
#    -k <value>   largest k-mer size to consider (default: 121)
#    -l <value>   smallest k-mer size to consider (default: 15)
#    -s <value>   interval between consecutive kmer sizes (default: 10)
#    -e <value>   k-mer sampling value (default: auto-detected to use ~200 MB memory/thread)
#    -t <value>   number of threads (default: number of cores minus one)
#    -o <prefix>  prefix of the output files (default: histograms)
