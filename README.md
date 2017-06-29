# Introduction
FMOE is an overlap-based error-correction algorithm using FM-index.

# Compile
To compile FMOE in your environment, type 

      1. ./autogen.sh 
      2. ./configure
      3. make -j 8

An executable program called FMOE will be found under the FMOE folder.

# Execution
FMOE tales fasta as input. If your reads are original Illumina PE reads in fastq format, you should run preprocess to discard quality and convert to fasta format.

	./FMOE preprocess --discard-quality -p 1 R1.fq R2.fq -o InputFile.fa
	./FMOE index -t 20 InputFile.fa
	./FMOE correct -t 20 -k 31 -K 15 InputFile.fa -o OutputFile.fa

InputFile.fa : Input file (fasta)
OutputFile.fa : the file of corrected read.


-k : the kmer size for identifying solid kmers during seeding stage.
-K : the kmer size when aligning compressed reads against the query during FM-index extension.
