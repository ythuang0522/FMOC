# Introduction
FMOE is an overlap-based error-correction algorithm using FM-index.

# Compile
FMOE utilizes zlib and Google Sparsehash libraries. After installing the two dependencies, compile FMOE by typing 

      1. ./autogen.sh 
      2. ./configure
      3. make -j 8

An executable program called FMOE will be found under the FMOE folder.

# Execution
FMOE tales fasta as input. If your reads are original Illumina PE reads in fastq format (e.g., R1.fq and R2.fq), you should run preprocess to discard quality and convert to fasta format.

	./FMOE preprocess --discard-quality -p 1 R1.fq R2.fq -o InputFile.fa
	./FMOE index -t 20 InputFile.fa
	./FMOE correct -t 20 -k 31 -K 15 InputFile.fa -o OutputFile.fa

InputFile.fa : Input file of raw reads

OutputFile.fa : Output file of corrected reads.

-k : the kmer size for identifying solid kmers during seeding stage.

-K : the kmer size when aligning compressed reads against the query during FM-index extension.
