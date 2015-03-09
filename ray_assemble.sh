#!/bin/bash
prefix=$1
bestk=`tail -n 1 log/kmergenie.$prefix.stdout | sed -e "s/best k: //"`

mpiexec -n 32 bin/Ray-v2.2.0/Ray \
-k $bestk \
-p ec_data/$prefix.paired.A.fastq ec_data/$prefix.paired.B.fastq \
-s ec_data/$prefix.unpaired.fastq \
-o assembly/ray.k.$bestk.$prefix

#NAME
#       Ray - assemble genomes in parallel using the message-passing interface
#
#SYNOPSIS
#       mpiexec -n 80 Ray -k 31 -p l1_1.fastq l1_2.fastq -p l2_1.fastq l2_2.fastq -o test
#
#       mpiexec -n 80 Ray Ray.conf # with commands in a file
#
#       mpiexec -n 10 Ray -mini-ranks-per-rank 7 Ray.conf # with mini-ranks
#
#DESCRIPTION:
#
#  The Ray genome assembler is built on top of the RayPlatform, a generic plugin-based
#  distributed and parallel compute engine that uses the message-passing interface
#  for passing messages.
#
#  Ray targets several applications:
#
#    - de novo genome assembly (with Ray vanilla)
#    - de novo meta-genome assembly (with Ray Méta)
#    - de novo transcriptome assembly (works, but not tested a lot)
#    - quantification of contig abundances
#    - quantification of microbiome consortia members (with Ray Communities)
#    - quantification of transcript expression
#    - taxonomy profiling of samples (with Ray Communities)
#    - gene ontology profiling of samples (with Ray Ontologies)
#
#       -help
#              Displays this help page.
#
#       -version
#              Displays Ray version and compilation options.
#
#  Run Ray in pure MPI mode
#
#    mpiexec -n 80 Ray ...
#
#  Run Ray with mini-ranks on 10 machines, 8 cores / machine (MPI and IEEE POSIX threads)
#
#    mpiexec -n 10 Ray -mini-ranks-per-rank 7 ...
#
#  Run Ray on one core only (still needs MPI)
#
#    Ray ...
#
#
#  Using a configuration file
#
#    Ray can be launched with
#    mpiexec -n 16 Ray Ray.conf
#    The configuration file can include comments (starting with #).
#
#  K-mer length
#
#       -k kmerLength
#              Selects the length of k-mers. The default value is 21. 
#              It must be odd because reverse-complement vertices are stored together.
#              The maximum length is defined at compilation by CONFIG_MAXKMERLENGTH
#              Larger k-mers utilise more memory.
#
#  Inputs
#
#       -p leftSequenceFile rightSequenceFile [averageOuterDistance standardDeviation]
#              Provides two files containing paired-end reads.
#              averageOuterDistance and standardDeviation are automatically computed if not provided.
#
#       -i interleavedSequenceFile [averageOuterDistance standardDeviation]
#              Provides one file containing interleaved paired-end reads.
#              averageOuterDistance and standardDeviation are automatically computed if not provided.
#
#       -s sequenceFile
#              Provides a file containing single-end reads.
#
#  Outputs
#
#       -o outputDirectory
#              Specifies the directory for outputted files. Default is RayOutput
#
#  Assembly options (defaults work well)
#
#       -disable-recycling
#              Disables read recycling during the assembly
#              reads will be set free in 3 cases:
#              1. the distance did not match for a pair
#              2. the read has not met its mate
#              3. the library population indicates a wrong placement
#              see Constrained traversal of repeats with paired sequences.
#              Sébastien Boisvert, Élénie Godzaridis, François Laviolette & Jacques Corbeil.
#              First Annual RECOMB Satellite Workshop on Massively Parallel Sequencing, March 26-27 2011, Vancouver, BC, Canada.
#
#       -debug-recycling
#              Debug the recycling events
#
#       -disable-scaffolder
#              Disables the scaffolder.
#
#       -minimum-contig-length minimumContigLength
#              Changes the minimum contig length, default is 100 nucleotides
#
#       -color-space
#              Runs in color-space
#              Needs csfasta files. Activated automatically if csfasta files are provided.
#
#       -use-maximum-seed-coverage maximumSeedCoverageDepth
#              Ignores any seed with a coverage depth above this threshold.
#              The default is 4294967295.
#
#       -use-minimum-seed-coverage minimumSeedCoverageDepth
#              Sets the minimum seed coverage depth.
#              Any path with a coverage depth lower than this will be discarded. The default is 0.
#
#  Distributed storage engine (all these values are for each MPI rank)
#
#       -bloom-filter-bits bits
#              Sets the number of bits for the Bloom filter
#              Default is auto bits (adaptive), 0 bits disables the Bloom filter.
#
#       -hash-table-buckets buckets
#              Sets the initial number of buckets. Must be a power of 2 !
#              Default value: 268435456
#
#       -hash-table-buckets-per-group buckets
#              Sets the number of buckets per group for sparse storage
#              Default value: 64, Must be between >=1 and <= 64
#
#       -hash-table-load-factor-threshold threshold
#              Sets the load factor threshold for real-time resizing
#              Default value: 0.75, must be >= 0.5 and < 1
#
#       -hash-table-verbosity
#              Activates verbosity for the distributed storage engine
#
#  Biological abundances
#
#       -search searchDirectory
#              Provides a directory containing fasta files to be searched in the de Bruijn graph.
#              Biological abundances will be written to RayOutput/BiologicalAbundances
#              See Documentation/BiologicalAbundances.txt
#
#       -one-color-per-file
#              Sets one color per file instead of one per sequence.
#              By default, each sequence in each file has a different color.
#              For files with large numbers of sequences, using one single color per file may be more efficient.
#
#  Taxonomic profiling with colored de Bruijn graphs
#
#       -with-taxonomy Genome-to-Taxon.tsv TreeOfLife-Edges.tsv Taxon-Names.tsv
#              Provides a taxonomy.
#              Computes and writes detailed taxonomic profiles.
#              See Documentation/Taxonomy.txt for details.
#
#       -gene-ontology OntologyTerms.txt  Annotations.txt
#              Provides an ontology and annotations.
#              OntologyTerms.txt is fetched from http://geneontology.org
#              Annotations.txt is a 2-column file (EMBL_CDS handle	&	gene ontology identifier)
#              See Documentation/GeneOntology.txt
#  Other outputs
#
#       -enable-neighbourhoods
#              Computes contig neighborhoods in the de Bruijn graph
#              Output file: RayOutput/NeighbourhoodRelations.txt
#
#       -amos
#              Writes the AMOS file called RayOutput/AMOS.afg
#              An AMOS file contains read positions on contigs.
#              Can be opened with software with graphical user interface.
#
#       -write-kmers
#              Writes k-mer graph to RayOutput/kmers.txt
#              The resulting file is not utilised by Ray.
#              The resulting file is very large.
#
#       -write-read-markers
#              Writes read markers to disk.
#
#       -write-seeds
#              Writes seed DNA sequences to RayOutput/Rank<rank>.RaySeeds.fasta
#
#       -write-extensions
#              Writes extension DNA sequences to RayOutput/Rank<rank>.RayExtensions.fasta
#
#       -write-contig-paths
#              Writes contig paths with coverage values
#              to RayOutput/Rank<rank>.RayContigPaths.txt
#
#       -write-marker-summary
#              Writes marker statistics.
#
#  Memory usage
#
#       -show-memory-usage
#              Shows memory usage. Data is fetched from /proc on GNU/Linux
#              Needs __linux__
#
#       -show-memory-allocations
#              Shows memory allocation events
#
#  Algorithm verbosity
#
#       -show-extension-choice
#              Shows the choice made (with other choices) during the extension.
#
#       -show-ending-context
#              Shows the ending context of each extension.
#              Shows the children of the vertex where extension was too difficult.
#
#       -show-distance-summary
#              Shows summary of outer distances used for an extension path.
#
#       -show-consensus
#              Shows the consensus when a choice is done.
#
#  Checkpointing
#
#       -write-checkpoints checkpointDirectory
#              Write checkpoint files
#
#       -read-checkpoints checkpointDirectory
#              Read checkpoint files
#
#       -read-write-checkpoints checkpointDirectory
#              Read and write checkpoint files
#
#  Message routing for large number of cores
#
#       -route-messages
#              Enables the Ray message router. Disabled by default.
#              Messages will be routed accordingly so that any rank can communicate directly with only a few others.
#              Without -route-messages, any rank can communicate directly with any other rank.
#              Files generated: Routing/Connections.txt, Routing/Routes.txt and Routing/RelayEvents.txt
#              and Routing/Summary.txt
#
#       -connection-type type
#              Sets the connection type for routes.
#              Accepted values are debruijn, hypercube, polytope, group, random, kautz and complete. Default is debruijn.
#               torus: a k-ary n-cube, radix: k, dimension: n, degree: 2*dimension, vertices: radix^dimension
#               polytope: a convex regular polytope, alphabet is {0,1,...,B-1} and the vertices is a power of B
#               hypercube: a hypercube, alphabet is {0,1} and the vertices is a power of 2
#               debruijn: a full de Bruijn graph a given alphabet and diameter
#               kautz: a full de Kautz graph, which is a subgraph of a de Bruijn graph
#               group: silly model where one representative per group can communicate with outsiders
#               random: Erdős–Rényi model
#               complete: a full graph with all the possible connections
#              With the type debruijn, the number of ranks must be a power of something.
#              Examples: 256 = 16*16, 512=8*8*8, 49=7*7, and so on.
#              Otherwise, don't use debruijn routing but use another one
#              With the type kautz, the number of ranks n must be n=(k+1)*k^(d-1) for some k and d
#
#       -routing-graph-degree degree
#              Specifies the outgoing degree for the routing graph.
#              See Documentation/Routing.txt
#
#  Hardware testing
#
#       -test-network-only
#              Tests the network and returns.
#
#       -write-network-test-raw-data
#              Writes one additional file per rank detailing the network test.
#
#       -exchanges NumberOfExchanges
#              Sets the number of exchanges
#
#       -disable-network-test
#              Skips the network test.
#
#  Debugging
#
#       -verify-message-integrity
#              Checks message data reliability for any non-empty message.
#              add '-D CONFIG_SSE_4_2' in the Makefile to use hardware instruction (SSE 4.2)
#
#       -run-profiler
#              Runs the profiler as the code runs. By default, only show granularity warnings.
#              Running the profiler increases running times.
#
#       -with-profiler-details
#              Shows number of messages sent and received in each methods during in each time slices (epochs). Needs -run-profiler.
#
#       -show-communication-events
#              Shows all messages sent and received.
#
#       -show-read-placement
#              Shows read placement in the graph during the extension.
#
#       -debug-bubbles
#              Debugs bubble code.
#              Bubbles can be due to heterozygous sites or sequencing errors or other (unknown) events
#
#       -debug-seeds
#              Debugs seed code.
#              Seeds are paths in the graph that are likely unique.
#
#       -debug-fusions
#              Debugs fusion code.
#
#       -debug-scaffolder
#              Debug the scaffolder.
#
#
#FILES
#
#  Input files
#
#     Note: file format is determined with file extension.
#
#     .fasta
#     .fasta.gz (needs HAVE_LIBZ=y at compilation)
#     .fasta.bz2 (needs HAVE_LIBBZ2=y at compilation)
#     .fastq
#     .fastq.gz (needs HAVE_LIBZ=y at compilation)
#     .fastq.bz2 (needs HAVE_LIBBZ2=y at compilation)
#     export.txt
#     .sff (paired reads must be extracted manually)
#     .csfasta (color-space reads)
#
#  Outputted files
#
#  Scaffolds
#
#     RayOutput/Scaffolds.fasta
#     	The scaffold sequences in FASTA format
#     RayOutput/ScaffoldComponents.txt
#     	The components of each scaffold
#     RayOutput/ScaffoldLengths.txt
#     	The length of each scaffold
#     RayOutput/ScaffoldLinks.txt
#     	Scaffold links
#
#  Contigs
#
#     RayOutput/Contigs.fasta
#     	Contiguous sequences in FASTA format
#     RayOutput/ContigLengths.txt
#     	The lengths of contiguous sequences
#
#  Summary
#
#     RayOutput/OutputNumbers.txt
#     	Overall numbers for the assembly
#
#  de Bruijn graph
#
#     RayOutput/CoverageDistribution.txt
#     	The distribution of coverage values
#     RayOutput/CoverageDistributionAnalysis.txt
#     	Analysis of the coverage distribution
#     RayOutput/degreeDistribution.txt
#     	Distribution of ingoing and outgoing degrees
#     RayOutput/kmers.txt
#     	k-mer graph, required option: -write-kmers
#         The resulting file is not utilised by Ray.
#         The resulting file is very large.
#
#  Assembly steps
#
#     RayOutput/SeedLengthDistribution.txt
#         Distribution of seed length
#     RayOutput/Rank<rank>.OptimalReadMarkers.txt
#         Read markers.
#     RayOutput/Rank<rank>.RaySeeds.fasta
#         Seed DNA sequences, required option: -write-seeds
#     RayOutput/Rank<rank>.RayExtensions.fasta
#         Extension DNA sequences, required option: -write-extensions
#     RayOutput/Rank<rank>.RayContigPaths.txt
#         Contig paths with coverage values, required option: -write-contig-paths
#
#  Paired reads
#
#     RayOutput/LibraryStatistics.txt
#     	Estimation of outer distances for paired reads
#     RayOutput/LibraryData.xml
#         Frequencies for observed outer distances (insert size + read lengths)
#
#  Partition
#
#     RayOutput/NumberOfSequences.txt
#         Number of reads in each file
#     RayOutput/SequencePartition.txt
#     	Sequence partition
#
#  Ray software
#
#     RayOutput/RayVersion.txt
#     	The version of Ray
#     RayOutput/RayCommand.txt
#     	The exact same command provided 
#
#  AMOS
#
#     RayOutput/AMOS.afg
#     	Assembly representation in AMOS format, required option: -amos
#
#  Communication
#
#     RayOutput/NetworkTest.txt
#	    	Latencies in microseconds
#     RayOutput/Rank<rank>NetworkTestData.txt
#	    	Network test raw data
#
#DOCUMENTATION
#
#       - mpiexec -n 1 Ray -help|less (always up-to-date)
#       - This help page (always up-to-date)
#       - The directory Documentation/
#       - Manual (Portable Document Format): InstructionManual.tex (in Documentation)
#       - Mailing list archives: http://sourceforge.net/mailarchive/forum.php?forum_name=denovoassembler-users
#
#AUTHOR
#       Written by Sébastien Boisvert.
#
#REPORTING BUGS
#       Report bugs to denovoassembler-users@lists.sourceforge.net
#       Home page: <http://denovoassembler.sourceforge.net/>
#
#COPYRIGHT
#       This program is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, version 3 of the License.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You have received a copy of the GNU General Public License
#       along with this program (see LICENSE).
#
#Ray 2.2.0
