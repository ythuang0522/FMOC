//Construct in 2016/1/11 by Yu-Wen

#ifndef EXTENSION_H
#define EXTENSION_H

#include <stdio.h>
#include "BWT.h"
#include "BWTAlgorithms.h"
#include "FMExObject.h"
#include "FMExtendTree.h"
#include "FMExtendTreeRC.h"



// functions
namespace Extension
{
	bool isLowComplexity (std::string seq, float threshold=0.9);
	
	void getLRKmerInterval(std::string Query,int kmer_size,const BWT* pBWT, const BWT* pRBWT,
					std::vector<BWTInterval>& LTInterval,std::vector<BWTInterval>& RTInterval);
	
	void addQueryinOutInfo(std::string Query,std::vector<OutInfo>& Out_Info);
	
	void addQueryinKsub(std::string Query,int kmer_size,std::vector<Ksub_vct>& correct_ksub);
	
	void addStrInKsub(std::vector<OutInfo>& Out_Info,int kmer_size,std::vector<Ksub_vct>& correct_ksub);
	void addStrInKsub2(std::vector<OutInfo>& Out_Info,int kmer_size,std::vector<Out_test>& out_vct,std::string Query);
	void addStrInKsub3(std::vector<OutInfo>& Out_Info,int kmer_size,std::vector<Ksub_vct>& correct_ksub,std::string Query);
	
	Solid_error getSolidRegion(std::string Query,int kmer_size,int thrshold,const BWT* pBWT);
	Solid_error highError_getSolidRegion(std::string Query,int kmer_size,const BWT* pBWT);
	std::string getSolidRegion_v2(std::string Query,int Seed_size,const BWT* pBWT);
	
	bool ExtensionRead(std::string Query,int kmer_size,int check_kmer_size,int solid_idx,const BWT* pBWT, const BWT* pRBWT,
					std::vector<OutInfo>& Out_Info,std::vector<BWTInterval>& LTInterval,std::vector<BWTInterval>& RTInterval);
					
	std::string correct_right(Solid_error& solid_info,std::vector<Ksub_vct>& correct_ksub,int k_diff,std::string last_kmer);
	std::string correct_left(Solid_error& solid_info,std::vector<Ksub_vct>& correct_ksub,std::string last_kmer);
	
	std::string Newcorrect_right(std::string Query,Solid_error& solid_info,std::vector<Ksub_vct>& correct_ksub,int k_diff,std::string last_kmer,std::string FirstKmer,int Solid_kmer_size,const BWT* pBWT);
	std::string Newcorrect_left(std::string Query,Solid_error& solid_info,std::vector<Ksub_vct>& correct_ksub,std::string last_kmer,std::string FirstKmer,int Solid_kmer_size,const BWT* pBWT);
	
	std::string TrimReads(std::string consensus,int kmer_size,const BWT* pBWT);
	std::string NewTrimReads(std::string consensus,int kmer_size,const BWT* pBWT);
	void printfKFQ(std::string Query,int kmer_size,const BWT* pBWT);
					
};

#endif
