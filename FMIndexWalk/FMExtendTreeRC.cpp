//FMExtendTreeRC.cpp
//FM-index extenson tree
///----------------------------------------------
#include "FMExtendTreeRC.h"
#include "BWTAlgorithms.h"
#include<string.h>


//std::vector<Mstr>& out_str,
//							   std::vector<BWTInterval>& L_TerminatedIntervals,
//							   std::vector<BWTInterval>& R_TerminatedIntervals,
//m_out_str(out_str)

// Class: FMExtendTreeRC
FMExtendTreeRC::FMExtendTreeRC(const std::string& query,
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
							   std::vector<BWTInterval>& L_TerminatedIntervals,//向左向右的答案都從左邊開始切
							   std::vector<BWTInterval>& R_TerminatedIntervals,
							   size_t maxLeaves,
							   size_t seedDist):
                               m_Query(query),m_Kmer(kmer), m_initSeedoffset(initSeedoffset),m_minOverlap(minOverlap),  
                               m_minIdentity(minIdentity),m_maxIndelSize(15),
							   m_pBWT(pBWT), m_pRBWT(pRBWT),m_Out_info(Out_Info),
							   m_L_TerminatedIntervals(L_TerminatedIntervals),m_R_TerminatedIntervals(R_TerminatedIntervals),
                               m_maxLeaves(maxLeaves), m_seedSize(check_kmer_size), m_seedDist(seedDist)
{
	// create one root node
	FMExOverlapNode* m_pRootNode = new FMExOverlapNode(&m_Kmer, NULL);
	m_pRootNode->currIntervalPair = bip;
	//m_pRootNode->rvc_currIntervalPair = rvc_bip;
	//size_t var_k = m_Kmer.length()-m_seedSize;
	
	m_pRootNode->SetInitial();   //store initial str of root
	//L_currentLength = m_Query.length()-m_initSeedoffset;//m_seedSize + m_initSeedoffset;
	//R_currentLength = m_initSeedoffset+m_Kmer.length();
	
	L_currentLength = m_initSeedoffset+m_Kmer.length();//m_seedSize + m_initSeedoffset;
	R_currentLength = m_Query.length()-m_initSeedoffset;
	
	//maybe is not use the variances below
	// ----m_seedSize+m_initSeedoffset;----
	//maybe is not use the variances upside
	
	//m_pRootNode->L_lastSeedIdx = m_Query.length()- m_seedSize- m_initSeedoffset;
	//m_pRootNode->R_lastSeedIdx = m_initSeedoffset+var_k;
	//m_pRootNode->L_nextSeedIdx = m_pRootNode->L_lastSeedIdx + 1;
	//m_pRootNode->R_nextSeedIdx = m_pRootNode->R_lastSeedIdx + 1;
	
	// initial set to left-direct extend
	//m_pRootNode->lastSeedIdx = m_pRootNode->initSeedIdx = m_initSeedoffset;
	//m_pRootNode->nextSeedIdx=m_pRootNode->lastSeedIdx-1;
	
	m_pRootNode->lastSeedIdx = m_pRootNode->initSeedIdx = m_initSeedoffset+(m_Kmer.length()-m_seedSize);
	m_pRootNode->nextSeedIdx=m_pRootNode->lastSeedIdx+1;
	
	//m_pRootNode->real_nowidx =(int)(m_pRootNode->lastSeedIdx+1);
	m_pRootNode->numOfMatchSeeds = 0;
			
	// push new node into roots and leaves vector
	m_RootNodes.push_back(m_pRootNode);
	m_leaves.push_back(m_pRootNode);

	
}

//
FMExtendTreeRC::~FMExtendTreeRC()
{
	for (std::list<FMExOverlapNode*>::iterator it = m_RootNodes.begin(); it != m_RootNodes.end(); ++it)
		delete *it;
	
	m_RootNodes.clear();
    // Recursively destroy the tree
    // delete m_pRootNode;
}

int FMExtendTreeRC::extendOneBase(int direction)
{
	if(direction==0)
	{
		//Overlap extension via FM-index walk
		if(!m_leaves.empty() && m_leaves.size() <= m_maxLeaves )
		{
			// ACGT-extend the leaf nodes via updating existing SA interval
			extendLeaves(0);	
			// Remove leaves without seed support within m_maxIndelSize
			//PrunedBySeedSupport(0);
			//see if any interval reach left end $ with sufficient overlap
			isTerminated_L();	
		}
				
		//Did not reach the terminal kmer
		if(m_leaves.empty())
			return -1;	//high error
		else if(m_leaves.size() > m_maxLeaves)
			return -2;	//too much repeats
		else
			return 1;
	}
	else if(direction==1)
	{
		//Overlap extension via FM-index walk
		if(!m_leaves.empty() && m_leaves.size() <= m_maxLeaves )
		{
			// ACGT-extend the leaf nodes via updating existing SA interval
			extendLeaves(1);
			// Remove leaves without seed support within m_maxIndelSize
			//PrunedBySeedSupport(1);
			//see if any interval reach left end $ with sufficient overlap
			isTerminated_R();	
		}
				
		//Did not reach the terminal kmer
		if(m_leaves.empty())
			return -1;	//high error
		else if(m_leaves.size() > m_maxLeaves)
			return -2;	//too much repeats
		else
			return 1;
	}
	return 1;
}

void FMExtendTreeRC::extendLeaves(int direction)
{
	if(direction==0)
	{
		FMEONodePtrList newLeaves;
		newLeaves.clear();
		//attempt to extend one base for each leave
		attempToExtend_L(newLeaves);

		L_currentLength++;
		m_leaves.clear();
		m_leaves = newLeaves;
	}
	else
	{
		FMEONodePtrList newLeaves;
		newLeaves.clear();
		//attempt to extend one base for each leave
		attempToExtend_R(newLeaves);
			
		R_currentLength++;
		m_leaves.clear();
		m_leaves = newLeaves;
	}
}

void FMExtendTreeRC::attempToExtend_L(FMEONodePtrList &newLeaves)
{
    for(FMEONodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
		//std::string out_s=reverse((*iter)->getLeftExtend())+m_Kmer;
		//printf("left_lable=\t%s\n",out_s.c_str());
		//add the condition when the left-extend is out of Query
		
		//if((int)(*iter)->nextSeedIdx <0)
		if( (int)(*iter)->nextSeedIdx > ((int)m_R_TerminatedIntervals.size()-1) )
		{
			//printf("<0--The left-lable add to left_temp>\n%s\n",out_s.c_str());
			Left_tmp_leaves.push_back((*iter));
		}
		else
		{
			std::vector< std::pair<std::string, BWTIntervalPair> > extensions;
			extensions = getLeftFMIndexExtensions(*iter);
			
			if(extensions.size() == 1)
			{
				// Single extension, do not branch
				//(*iter)->extend_L(extensions.front().first);
				(*iter)->extend(extensions.front().first,0);
				(*iter)->currIntervalPair=extensions.front().second;
				(*iter)->nextSeedIdx++;
				// currOverlapLen/queryOverlapLen always increase wrt each extension
				// in order to know the approximate real-time matched length for terminal/containment processing
				//(*iter)->currOverlapLen++;
				//(*iter)->queryOverlapLen++;
				//(*iter)->node_currlength++;
				newLeaves.push_back(*iter);
			}
			else if(extensions.size() > 1)
			{
				// Branch for each child
				for(size_t i = 0; i < extensions.size(); ++i)
				{
					FMExOverlapNode* pChildNode = (*iter)->createChild(extensions[i].first,0);
					pChildNode->currIntervalPair = extensions[i].second;
					pChildNode->nextSeedIdx++;
					//pChildNode->currOverlapLen++;
					//pChildNode->queryOverlapLen++;
					//pChildNode->node_currlength++;
					
					newLeaves.push_back(pChildNode);
				}
			}
		}
	}	
}

void FMExtendTreeRC::attempToExtend_R(FMEONodePtrList &newLeaves)
{
	//printf("The_Right_length\t%d\n",(int)R_currentLength);
    for(FMEONodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
		//std::string out_s=(*iter)->getRightExtend();
		//printf("rigth_lable=\t%s\n",out_s.c_str());
		//add the condition when the right-extend is out of Query
		
		//if( (int)(*iter)->nextSeedIdx > ((int)m_R_TerminatedIntervals.size()-1) )
		if((int)(*iter)->nextSeedIdx <0)
		{
			OutToCorrect((*iter));
			/*
			int left_gap=(int)m_initSeedoffset-(int)(*iter)->getLeftExtend().length();
			std::string l_gap_str="";
			for(int t=0;t<left_gap;t++)
			{
				l_gap_str.append("=");
			}
			
			std::string out_extend_str=l_gap_str+reverse((*iter)->getLeftExtend())+"~"+m_Kmer+"~"+(*iter)->getRightExtend();
			//printf("Ex_str=\t%s-%d\n",out_extend_str.c_str(),(int)(*iter)->nextSeedIdx);
			printf("Ex_str=\t%s-%d\t%d\nMis=%d\tdele=%d\tInsert=%d\ttempMis=%d\n",out_extend_str.c_str(),(int)(*iter)->nextSeedIdx,(int)(*iter)->currIntervalPair.interval[1].lower,(*iter)->mismatch_count,
			(*iter)->deletion_count,(*iter)->insertion_count,(*iter)->tempMismathCount);
			if((*iter)->insertion_count>0)
			{
				for( int i=0;i<(int)(*iter)->Insertion_vector.size();i++ )
				{
					printf("Idx=%d\tpreKmer=%s\tinsert_str=%s\n",(*iter)->Insertion_vector[i].match_idx,
					(*iter)->Insertion_vector[i].preKmer.c_str(),(*iter)->Insertion_vector[i].insertion_str.c_str());
				}
			}
			*/
		}
		else
		{
			std::vector< std::pair<std::string, BWTIntervalPair> > extensions;
			extensions = getRightFMIndexExtensions(*iter);
			
			if(extensions.size() == 1)
			{
				// Single extension, do not branch
				(*iter)->extend(extensions.front().first,1);
				(*iter)->currIntervalPair=extensions.front().second;
				(*iter)->nextSeedIdx--;
				// currOverlapLen/queryOverlapLen always increase wrt each extension
				// in order to know the approximate real-time matched length for terminal/containment processing
				//(*iter)->currOverlapLen++;
				//(*iter)->queryOverlapLen++;
				//(*iter)->node_currlength++;
				
				newLeaves.push_back(*iter);
			}
			else if(extensions.size() > 1)
			{
				// Branch for each child
				for(size_t i = 0; i < extensions.size(); ++i)
				{
					FMExOverlapNode* pChildNode = (*iter)->createChild(extensions[i].first,1);
					pChildNode->currIntervalPair = extensions[i].second;
					pChildNode->nextSeedIdx--;
					//pChildNode->currOverlapLen++;
					//pChildNode->queryOverlapLen++;
					//pChildNode->node_currlength++;
					
					newLeaves.push_back(pChildNode);
				}
			}
		}
    }	
}

std::vector<std::pair<std::string, BWTIntervalPair> > FMExtendTreeRC::getLeftFMIndexExtensions(FMExOverlapNode* pNode)
{
    std::vector<std::pair<std::string, BWTIntervalPair> > out;

    for(int i = 1; i < BWT_ALPHABET::size; ++i) //i=A,C,G,T
    {
        char b = BWT_ALPHABET::getChar(i);

        //update IntervalPair using extension b
        BWTIntervalPair probe=pNode->currIntervalPair;
        BWTAlgorithms::updateBothL(probe, b, m_pBWT);
		
		//char rvc_b = BWT_ALPHABET::getChar(5-i);
        //update reverse complement IntervalPair using extension b
        //BWTIntervalPair rvc_probe=pNode->rvc_currIntervalPair;
        //BWTAlgorithms::updateBothL(rvc_probe, rvc_b, m_pRBWT);
		//min freq at fwd and rvc bwt
        if(probe.isValid())
        {
            std::string tmp;
            tmp.append(1,b);
            BWTIntervalPair bip=probe;
            out.push_back(std::make_pair(tmp, bip));
        }
    }// end of ACGT

    return out;
}
std::vector<std::pair<std::string, BWTIntervalPair> > FMExtendTreeRC::getRightFMIndexExtensions(FMExOverlapNode* pNode)
{
    std::vector<std::pair<std::string, BWTIntervalPair> > out;

    for(int i = 1; i < BWT_ALPHABET::size; ++i) //i=A,C,G,T
    {
        char b = BWT_ALPHABET::getChar(i);

        //update IntervalPair using extension b
        BWTIntervalPair probe=pNode->currIntervalPair;
        BWTAlgorithms::updateBothR(probe, b, m_pRBWT);
		
		//char rvc_b = BWT_ALPHABET::getChar(5-i);
        //update reverse complement IntervalPair using extension b
       // BWTIntervalPair rvc_probe=pNode->rvc_currIntervalPair;
        //BWTAlgorithms::updateBothR(rvc_probe, rvc_b, m_pBWT);		
		//min freq at fwd and rvc bwt
        if( probe.isValid())
        {
            std::string tmp;
            tmp.append(1,b);
            BWTIntervalPair bip=probe;
            out.push_back(std::make_pair(tmp, bip));
        }
    }// end of ACGT
    return out;
}
/*
bool FMExtendTreeRC::isEnoughOverlap(FMExOverlapNode* pNode)
{
	int p_Lidx = (int)(m_Query.length()-m_seedSize-(pNode->L_lastSeedIdx));
	int p_Ridx = (int)( (pNode->R_lastSeedIdx)+m_seedSize-1 );
	int OlpLen = p_Ridx - p_Lidx;
	if( OlpLen >= (int)m_Query.length()/3 )
		return true;
	else
		return false;
}
*/
/*
void FMExtendTreeRC::OutResult(FMExOverlapNode* pNode)
{
	int flag=1;
	for(int i=0;i<(int)m_out_str.size();i++)
	{
		if( m_out_str[i].lower<= (int)pNode->currIntervalPair.interval[1].lower &&
		m_out_str[i].upper >= (int)pNode->currIntervalPair.interval[1].upper)
		{ 
			flag=0;
			break;
		}
	}
	if(flag==1)
	{
		std::string tmp_str=adjust_byIndel(pNode);
		if( tmp_str.length() >= m_minOverlap)
		{
			//std::string tmp_str=pNode->printLeft()+pNode->printRight();
			Mstr temp_mstr;
			temp_mstr.extend_str = tmp_str;
			temp_mstr.interval_size = (int)pNode->currIntervalPair.interval[1].size();
			temp_mstr.lower = (int)pNode->currIntervalPair.interval[1].lower;
			temp_mstr.upper = (int)pNode->currIntervalPair.interval[1].upper;
			temp_mstr.Left_start_idx = pNode->Left_start_idx;
			//temp_mstr.Right_end_idx = pNode->Right_end_idx=(int)(R_currentLength)-1;
			temp_mstr.Right_end_idx = pNode->Right_end_idx=(int)(R_currentLength)+pNode->de_count-pNode->in_count-1;
			temp_mstr.idl_n.insert=pNode->in_count;
			temp_mstr.idl_n.delet=pNode->de_count;
			temp_mstr.idl_n.mis=pNode->mis_count;
			temp_mstr.idl_n.tmpmis=pNode->tempMismatch;
			m_out_str.push_back(temp_mstr);
		}
	}
		
}
*/

// Check for leaves whose extension has terminated. If the leaf has
// terminated, the walked string and coverage is pushed to the result vector
bool FMExtendTreeRC::PrunedBySeedSupport(int directflag)
{
	FMEONodePtrList newLeaves;
	//int currSeedIdx;
	//directflag 0=Left 1=Right
	/*
	if(directflag==1)
		currSeedIdx = (R_currentLength-m_seedSize) / m_seedDist;
	else
		currSeedIdx = (L_currentLength-m_seedSize) / m_seedDist;
	*/
	//size_t indelOffset = (m_seedSize+m_maxIndelSize) / m_seedDist;
	size_t indelOffset = 5;
	//size_t smallSeedIdx = currSeedIdx <= indelOffset ? 0 : currSeedIdx - indelOffset;
	//size_t largeSeedIdx = (currSeedIdx+indelOffset) >= (m_L_TerminatedIntervals.size()-1)? 
	//						m_L_TerminatedIntervals.size()-1:currSeedIdx+indelOffset;

	// check range of last seed and find new seeds for each interval
    for(FMEONodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {	
		size_t currSeedIdx=(*iter)->nextSeedIdx;
		size_t largeSeedIdx;
		size_t smallSeedIdx;
		if(directflag==0)
		{
			smallSeedIdx = (int)(currSeedIdx+indelOffset) <= (int)(*iter)->lastSeedIdx ? (currSeedIdx+indelOffset) :(*iter)->lastSeedIdx;
			largeSeedIdx = (int)(currSeedIdx-indelOffset) >= 0 ? currSeedIdx-indelOffset:0;
		}
		else
		{
			smallSeedIdx = (int)(currSeedIdx-indelOffset) >= (int)(*iter)->lastSeedIdx ? (currSeedIdx-indelOffset) :(*iter)->lastSeedIdx;
			largeSeedIdx = (currSeedIdx+indelOffset) >= (m_R_TerminatedIntervals.size()-1)? m_L_TerminatedIntervals.size()-1:currSeedIdx+indelOffset;
		}
						
		bool isNewSeedFound;
		if(directflag==0)
			isNewSeedFound=isSupportedByNewSeed_L(*iter, smallSeedIdx,largeSeedIdx);
		else
			isNewSeedFound=isSupportedByNewSeed_R(*iter,smallSeedIdx, largeSeedIdx);

		
		bool indel_point=((*iter)->insertion_count+(*iter)->deletion_count + (*iter)->mismatch_count ) <(int)m_maxIndelSize;
		bool isOutofRange = (*iter)->tempMismathCount <= 10*(int)(indelOffset+m_seedSize);
		
		if(!(*iter)->isOverthreshold)
		{
			if( !( (*iter)->tempMismathCount <=(int)(indelOffset+m_seedSize) ))
			{
				int now_length = (int)(*iter)->getLeftExtend().length()+(int)(*iter)->getRightExtend().length();
				//printf("TempMismatch too much.\tNow Seeds=\t%d\tLength=\t%d\n",(int)(*iter)->numOfMatchSeeds,now_length);
				//printf("TempMuch-Seed_rate=\t%d/%d=%.2f\tlength=\t%d\nLeft=\t%s\nRight=\t%s\n",(int)(*iter)->numOfMatchSeeds,(int)(*iter)->getRightExtend().length(),((double)(*iter)->numOfMatchSeeds/(*iter)->getRightExtend().length()),now_length,(reverse((*iter)->getLeftExtend())+m_Kmer).c_str(),(*iter)->getRightExtend().c_str());
				(*iter)->seed_rate=(double)(*iter)->numOfMatchSeeds/(double)now_length;
				(*iter)->isOverthreshold=true;
			}
			if( !indel_point )
			{
				int now_length = (int)(*iter)->getLeftExtend().length()+(int)(*iter)->getRightExtend().length();
				//printf("ErrorMuch-Seed_rate=\t%d/%d=%.2f\tlength=\t%d\nLeft=\t%s\nRight=\t%s\n",(int)(*iter)->numOfMatchSeeds,(int)(*iter)->getRightExtend().length(),((double)(*iter)->numOfMatchSeeds/(*iter)->getRightExtend().length()),now_length,(reverse((*iter)->getLeftExtend())+m_Kmer).c_str(),(*iter)->getRightExtend().c_str());
				//printf("Indel/Mismatch too much.\tNow Seeds=\t%d\tLength=\t%d\n",(int)(*iter)->numOfMatchSeeds,now_length);
				(*iter)->seed_rate=(double)(*iter)->numOfMatchSeeds/(double)now_length;
				(*iter)->isOverthreshold=true;
			}
		}
		
		//if( (isNewSeedFound) || indel_point)
		if( (isOutofRange || isNewSeedFound) && indel_point)
		{
			newLeaves.push_back(*iter);
		}
		else
		{
		
			OutToCorrect((*iter));
			
		}
		
    }
	
	// update m_leaves with newLeaves;
	bool updated = false;
	if(m_leaves.size() != newLeaves.size())
	{
		m_leaves = newLeaves;
		updated=true;
	}

    return updated;
}

// Identify new seeds wrt currSeedIdx
bool FMExtendTreeRC::isSupportedByNewSeed_L(FMExOverlapNode* currNode,size_t smallSeedIdx, size_t largeSeedIdx)
{

		bool isNewSeedFound = false;
		
		BWTIntervalPair currIntervalPair = currNode->currIntervalPair;
		
		//int nowidx=currNode->L_nextSeedIdx;
		int IdxOfExtendBp=(int)currNode->nextSeedIdx;
		int start_idx=((int)smallSeedIdx);
		start_idx=(int)currNode->lastSeedIdx-1;
		for(int i=start_idx;i>=(int)largeSeedIdx;i--)
		{
			
			isNewSeedFound = currIntervalPair.interval[0].lower >= m_L_TerminatedIntervals.at(i).lower
						&& currIntervalPair.interval[0].upper <= m_L_TerminatedIntervals.at(i).upper;
			//isNewSeedFound=true;

			if(isNewSeedFound)
			{
				currNode->lastSeedIdx = i;
				currNode->nextSeedIdx = i-1;
				//currNode->numOfMatchSeeds++;
				
				if((int)currNode->lastSeedIdx<IdxOfExtendBp)//deletion
				{
					int NumOfDeletion=0;
					NumOfDeletion=IdxOfExtendBp-(int)currNode->lastSeedIdx;
					
					currNode->deletion_count+=NumOfDeletion;
					currNode->tempMismathCount=0;
					//printf("Deletion_start\n");
					currNode->adjust_byDeletion(0,(int)m_seedSize,NumOfDeletion);
					//printf("Deletion_end\n");
					//printf("Left_deletion\n");
				}
				if((int)currNode->lastSeedIdx==IdxOfExtendBp)//mismatch
				{
					if(currNode->tempMismathCount>0)
					{
						int NumOfMismatch=0;
						NumOfMismatch=(currNode->tempMismathCount+1-m_seedSize);
						
						currNode->mismatch_count+=NumOfMismatch;
						
						//indel mis_temp( (int)currNode->lastSeedIdx+(int)m_seedSize,NumOfMismatch);
						
						//indel mis_temp(m_Query.length()-(m_seedSize+(int)currNode->lastSeedIdx)+1,(currNode->tempMismatch-m_seedSize));
						//currNode->mismatch_vector.push_back(mis_temp);
						currNode->tempMismathCount=0;
					}
				}
				if((int)currNode->lastSeedIdx>IdxOfExtendBp)//Insertion
				{
					
					int NumOfInsert=0;
					NumOfInsert=(int)currNode->lastSeedIdx-IdxOfExtendBp;
					//if(NumOfInsert>5)
					if((int)currNode->getLeftExtend().length()-(int)m_seedSize-NumOfInsert<0)
					{
						currNode->tempMismathCount+=1;
					}
					else
					{
					currNode->insertion_count+=NumOfInsert;
					currNode->tempMismathCount=0;
					std::pair<std::string,std::string> preK_Istr;
					//printf("L-Insertion_start\n");
					preK_Istr=currNode->adjust_byInsertion(0,m_seedSize,NumOfInsert);
					//printf("L-Insertion_end\n");
					//printf("Left_Insert\tsize=\t%d\n",NumOfInsert);
					Insertion Insert_temp((int)currNode->lastSeedIdx,(int)currIntervalPair.interval[0].size(),preK_Istr.first,preK_Istr.second);
					currNode->Insertion_vector.push_back(Insert_temp);
					}
					//indel de_temp(m_Query.length()-(nowidx)-deletion_temp,(int)currNode->L_lastSeedIdx-nowidx);
					//currNode->deletion_vct.push_back(de_temp);
					
				}
				
				break;
			}
		}
		
		if(!isNewSeedFound)
		{
			currNode->tempMismathCount+=1;
			currNode->nextSeedIdx-=1;
		}

		return isNewSeedFound;
}

// Identify new seeds wrt currSeedIdx
bool FMExtendTreeRC::isSupportedByNewSeed_R(FMExOverlapNode* currNode,size_t smallSeedIdx, size_t largeSeedIdx)
{
bool isNewSeedFound = false;
		
		BWTIntervalPair currIntervalPair = currNode->currIntervalPair;
		
		//printf("Last=\t%d\tsmall=\t%d\n",(int)currNode->lastSeedIdx,(int)smallSeedIdx);
		//int nowidx=currNode->L_nextSeedIdx;
		int IdxOfExtendBp=(int)currNode->nextSeedIdx;
		int start_idx=((int)smallSeedIdx);
		for(int i=start_idx;i<=(int)largeSeedIdx;i++)
		{
			isNewSeedFound = currIntervalPair.interval[1].lower >= m_R_TerminatedIntervals.at(i).lower
						&& currIntervalPair.interval[1].upper <= m_R_TerminatedIntervals.at(i).upper;
			
			if(isNewSeedFound)
			{
				
				currNode->lastSeedIdx = i;
				currNode->nextSeedIdx = i+1;
				currNode->numOfMatchSeeds++;
				
				if((int)currNode->lastSeedIdx>IdxOfExtendBp)//deletion
				{
					int NumOfDeletion=0;
					NumOfDeletion=(int)currNode->lastSeedIdx-IdxOfExtendBp;
					
					currNode->deletion_count+=NumOfDeletion;
					currNode->tempMismathCount=0;
					//printf("Deletion_start\n");
					currNode->adjust_byDeletion(1,(int)m_seedSize,NumOfDeletion);
					//printf("Deletion_end\n");
					//printf("Left_deletion\n");
				}
				if((int)currNode->lastSeedIdx==IdxOfExtendBp)//mismatch
				{
					if(currNode->tempMismathCount>0)
					{
						int NumOfMismatch=0;
						NumOfMismatch=(currNode->tempMismathCount+1-m_seedSize);
						
						currNode->mismatch_count+=NumOfMismatch;
						
						//indel mis_temp( (int)currNode->lastSeedIdx+(int)m_seedSize,NumOfMismatch);
						
						//indel mis_temp(m_Query.length()-(m_seedSize+(int)currNode->lastSeedIdx)+1,(currNode->tempMismatch-m_seedSize));
						//currNode->mismatch_vector.push_back(mis_temp);
						currNode->tempMismathCount=0;
					}
				}
				if((int)currNode->lastSeedIdx<IdxOfExtendBp)//Insertion
				{
					//int R_length=(int)currNode->getRightExtend().length();
					int NumOfInsert=0;
					NumOfInsert=IdxOfExtendBp-(int)currNode->lastSeedIdx;
					//if(NumOfInsert>5)
					if((int)currNode->getRightExtend().length()-(int)m_seedSize-NumOfInsert<0)
					{
					
						//printf("Right_Insert\tsize=\t%d\tlength=\t%d\n",NumOfInsert,R_length);
						//printf("Next=\t%d\tLast=\t%d\tnow=\t%d\n",(int)currNode->nextSeedIdx,(int)currNode->lastSeedIdx,IdxOfExtendBp);
						currNode->tempMismathCount+=1;
					}
					else
					{
					//printf("Last=\t%d\tNow=\t%d\n",(int)currNode->lastSeedIdx,IdxOfExtendBp);
					currNode->insertion_count+=NumOfInsert;
					currNode->tempMismathCount=0;
					std::pair<std::string,std::string> preK_Istr;
					//printf("R-Insertion_start\n");
					preK_Istr=currNode->adjust_byInsertion(1,m_seedSize,NumOfInsert);
					//printf("R-Insertion_end\n");
					//printf("Right_Insert\tsize=\t%d\tlength=\t%d\n",NumOfInsert,R_length);
					Insertion Insert_temp((int)currNode->lastSeedIdx-m_seedSize,(int)currIntervalPair.interval[1].size(),preK_Istr.first,preK_Istr.second);
					currNode->Insertion_vector.push_back(Insert_temp);
					}
					//indel de_temp(m_Query.length()-(nowidx)-deletion_temp,(int)currNode->L_lastSeedIdx-nowidx);
					//currNode->deletion_vct.push_back(de_temp);
					
				}
				
				break;
			}
		}
		
		if(!isNewSeedFound)
		{
			currNode->tempMismathCount+=1;
			currNode->nextSeedIdx+=1;
		}

		return isNewSeedFound;
}

bool FMExtendTreeRC::isTerminated_L()
{
	bool found=false;
	FMEONodePtrList newLeaves;
	
    for(FMEONodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
		
		BWTIntervalPair probe = (*iter)->currIntervalPair;
		BWTAlgorithms::updateBothL(probe, '$', m_pBWT);
		
		// The SA interval reach $
		if( probe.isValid() )
		{
			//if((int)(*iter)->nextSeedIdx>=0)
			if((int)(*iter)->nextSeedIdx<=(int)(m_R_TerminatedIntervals.size()-m_seedSize))
			{
				newLeaves.push_back((*iter));
			}
			
			//std::string out_s=reverse((*iter)->getLeftExtend())+m_Kmer;
			//printf("$--The left-lable add to left_temp>\n%s\n",out_s.c_str());
			
			
			
			FMExOverlapNode* pChildNode = (*iter)->createChild("",0);
			pChildNode->currIntervalPair = probe;
			pChildNode->tempMismathCount=0;
			Left_tmp_leaves.push_back(pChildNode);
		}
		else{
			newLeaves.push_back((*iter));
		}
		
	}
	m_leaves.clear();
	m_leaves = newLeaves;
	return found;

}

bool FMExtendTreeRC::isTerminated_R()
{
	bool found=false;
	FMEONodePtrList newLeaves;
	
    for(FMEONodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
		
		BWTIntervalPair probe = (*iter)->currIntervalPair;
		BWTAlgorithms::updateBothR(probe, '$', m_pRBWT);
		
		// The SA interval reach $
		if( probe.isValid() )
		{
			//if((int)(*iter)->nextSeedIdx<=(int)(m_R_TerminatedIntervals.size()-m_seedSize))
			if((int)(*iter)->nextSeedIdx>=0)
			{
				newLeaves.push_back((*iter));
			}
			
			
			OutToCorrect((*iter));
			/*
			int left_gap=(int)m_initSeedoffset-(int)(*iter)->getLeftExtend().length();
			std::string l_gap_str="";
			for(int t=0;t<left_gap;t++)
			{
				l_gap_str.append("=");
			}
			std::string out_extend_str=l_gap_str+reverse((*iter)->getLeftExtend())+"~"+m_Kmer+"~"+(*iter)->getRightExtend();
			printf("Ex_str=\t%s$\t%d\nMis=%d\tdele=%d\tInsert=%d\ttempMis=%d\n",out_extend_str.c_str(),(int)(*iter)->currIntervalPair.interval[1].lower,(*iter)->mismatch_count,
			(*iter)->deletion_count,(*iter)->insertion_count,(*iter)->tempMismathCount);
			if((*iter)->insertion_count>0)
			{
				for( int i=0;i<(int)(*iter)->Insertion_vector.size();i++ )
				{
					printf("Idx=%d\tpreKmer=%s\tinsert_str=%s\n",(*iter)->Insertion_vector[i].match_idx,
					(*iter)->Insertion_vector[i].preKmer.c_str(),(*iter)->Insertion_vector[i].insertion_str.c_str());
				}
			}
			*/
			
			
			
		}
		else{
			newLeaves.push_back((*iter));
		}
		
	}
	m_leaves.clear();
	m_leaves = newLeaves;
	return found;

}
void FMExtendTreeRC::OutToCorrect(FMExOverlapNode* currNode)
{
	
	if(!currNode->isOverthreshold)
	{
		int now_length = (int)currNode->getLeftExtend().length()+(int)currNode->getRightExtend().length();
		//printf("TempMismatch too much.\tNow Seeds=\t%d\tLength=\t%d\n",(int)(*iter)->numOfMatchSeeds,now_length);
		
		currNode->seed_rate=(double)currNode->numOfMatchSeeds/(double)now_length;
		currNode->isOverthreshold=true;
	}
	
	std::string out_string=reverse(currNode->getLeftExtend())+m_Kmer+currNode->getRightExtend();
	out_string=reverseComplement(out_string);
	//printf(">RC_str\n%s\n",out_string.c_str());
	//int start_point=(int)m_initSeedoffset-(int)currNode->getLeftExtend().length();
	int start_point=(int)m_initSeedoffset-(int)currNode->getRightExtend().length();
	
	
	int temp_size=(int)currNode->currIntervalPair.interval[0].size();
	
	OutInfo temp_outinfo(out_string,temp_size,start_point,currNode->seed_rate);
	m_Out_info.push_back(temp_outinfo);
	/*
	for(int i=0;i<=(int)(out_string.length()-m_seedSize);i++)
	{
		bool isExist=false;
		//printf("temp_kmer_start\n");
		std::string temp_kmer=out_string.substr(i,m_seedSize);
		//printf("temp_kmer=%s\n",temp_kmer.c_str());
		int temp_size=(int)currNode->currIntervalPair.interval[0].size();
		for(int t=0;t<(int)m_correct_ksub[i+start_point].size();t++)
		{
			//printf("T_size=\t%d\n",t);
			if(m_correct_ksub[i+start_point][t].kmer.compare(temp_kmer)==0)
			{
				isExist=true;
				m_correct_ksub[i+start_point][t].countOfkmer+=temp_size;
				break;
			}
		}
		
		if(!isExist)
		{
			Ksubstr temp_ksub(temp_kmer,temp_size);
			//Ksub_vct temp_ksub_vct;
			m_correct_ksub[i+start_point].push_back(temp_ksub);
		}
		
	}
	*/
}
void FMExtendTreeRC::Left2Right()
{
	LDtoRD();
}
void FMExtendTreeRC::LDtoRD()
{
	m_leaves.clear();
	m_leaves = Left_tmp_leaves;
	//printf("The_left_temp_size=\t%d\n",(int)m_leaves.size());
	for(FMEONodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
		//printf("Left=%s\tInterval=\t%d\t%d\tInterval_size=\t%d\n",(*iter)->getLeftExtend().c_str(),(int)(*iter)->currIntervalPair.interval[0].lower,(int)(*iter)->currIntervalPair.interval[1].lower,(int)(*iter)->currIntervalPair.interval[0].size());
		
		(*iter)->lastSeedIdx=m_initSeedoffset;
		(*iter)->nextSeedIdx=(*iter)->lastSeedIdx-1;
	}
	
	/*
	for(FMEONodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
	{
		printf("L=\t%s\n",(*iter)->printLeft().c_str());
	}
	*/
}

