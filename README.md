# Introduction
FMOC is an overlap-based error-correction algorithm using FM-index. The algoritm identifies overlapping reads w.r.t a query read using the seed-and-extend alignment aginst compressed reads. Errors are then corrected using major k-mers from overlapping reads.

# Compile
FMOC utilizes zlib and Google Sparsehash libraries. After installing the two dependencies, compile FMOC by typing 

      1. ./autogen.sh 
      2. ./configure
      3. make -j 8

An executable program called FMOC will be found under the FMOC folder.

# Execution
FMOC tales fasta as input. If your reads are original Illumina PE reads in fastq format (e.g., R1.fq and R2.fq), you should run preprocess to discard quality and convert to fasta format.

	./FMOC preprocess --discard-quality -p 1 R1.fq R2.fq -o InputFile.fa
	./FMOC index -t 20 InputFile.fa
	./FMOC correct -t 20 -k 31 -K 15 InputFile.fa -o OutputFile.fa

InputFile.fa : Input file of raw reads

OutputFile.fa : Output file of corrected reads.

-k : the kmer size for identifying solid kmers during seeding stage.

-K : the kmer size when aligning compressed reads against the query during FM-index extension.

-t : number of threads.
