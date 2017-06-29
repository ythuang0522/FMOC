//FMExtendTreeRC.h
//FM-index extension tree
//----------------------------------------------

#ifndef FMExtendTreeRCRC_H
#define FMExtendTreeRCRC_H

#include <list>
#include "BWT.h"
#include "BWTAlgorithms.h"
#include "FMExtendNode.h"
#include "OverlapBlock.h"

	//std::vector<Mstr>& out_str,
		//			   std::vector<BWTInterval>& L_TerminatedIntervals,
			//		   std::vector<BWTInterval>& R_TerminatedIntervals,
class FMExtendTreeRC
{
    public:
        FMExtendTreeRC(const std::string& query,
					   const std::string& kmer,
					   size_t check_kmer_size,
					   BWTIntervalPair bip,
					   //BWTIntervalPair rvc_bip,
					   size_t initSeedoffset,
                       size_t minOverlap,
					   size_t minIdentity,
					   const BWT* pBWT, 
					   const BWT* pRBWT,
					   std::vector<OutInfo>& Out_Info,					   
					   std::vector<BWTInterval>& L_TerminatedIntervals,
					   std::vector<BWTInterval>& R_TerminatedIntervals,
                       size_t maxLeaves=128,
					   size_t seedDist=1);
		
        ~FMExtendTreeRC();
		
		//int extendOverlapOneBase_L();
		//int extendOverlapOneBase_R();
		int extendOneBase(int direction);
		
		inline size_t getCurrentLength_L(){return L_currentLength;};
		inline size_t getCurrentLength_R(){return R_currentLength;};
		
		// return emptiness of leaves
		inline bool isEmpty(){return m_leaves.empty();};
		
		void Left2Right();
		

    private:
        // Functions
        void LDtoRD();
		//bool isEnoughOverlap(FMExOverlapNode* pNode);
		//void OutResult(FMExOverlapNode* pNode);
		
        //void extendLeaves_L();
		//void extendLeaves_R();
		void extendLeaves(int direction);
		
		bool PrunedBySeedSupport(int direction);
		
		void attempToExtend_L(FMEONodePtrList &newLeaves);
		void attempToExtend_R(FMEONodePtrList &newLeaves);
		
		std::vector<std::pair<std::string, BWTIntervalPair> > getLeftFMIndexExtensions(FMExOverlapNode* pNode);
		std::vector<std::pair<std::string, BWTIntervalPair> > getRightFMIndexExtensions(FMExOverlapNode* pNode);
		
		
		
		bool isSupportedByNewSeed_L(FMExOverlapNode* currNode, size_t smallSeedIdx,size_t largeSeedIdx);
		bool isSupportedByNewSeed_R(FMExOverlapNode* currNode, size_t smallSeedIdx,size_t largeSeedIdx);
		
		bool isTerminated_L();
		bool isTerminated_R();
		
		void OutToCorrect(FMExOverlapNode* currNode);
		
		//待調整
		//std::string adjust_byIndel(FMExOverlapNode* currNode);
		
		
        //
        // Data
        //
        const std::string m_Query;
		const std::string m_Kmer;
        size_t m_initSeedoffset;
		size_t m_minOverlap;
		size_t m_minIdentity;
		size_t m_maxIndelSize;
		//double m_errorRate;
		const BWT* m_pBWT;
		const BWT* m_pRBWT;
		//std::vector<Mstr>& 	m_out_str;
		//std::vector<Ksub_vct>& m_correct_ksub;
		std::vector<OutInfo>& m_Out_info;
		std::vector<BWTInterval>& m_L_TerminatedIntervals;  
		std::vector<BWTInterval>& m_R_TerminatedIntervals; 
		// Optional parameters
        size_t m_maxLeaves;
		size_t m_seedSize;
		size_t m_seedDist;

        FMEONodePtrList m_leaves;
		FMEONodePtrList Left_tmp_leaves;

        // FMExOverlapNode* m_pRootNode;
        FMEONodePtrList m_RootNodes;

        //size_t m_currentLength; maybe not neccessary
		size_t L_currentLength;
		size_t R_currentLength;

        //std::vector<BWTInterval> m_TerminatedIntervals;
		 		
};

#endif
