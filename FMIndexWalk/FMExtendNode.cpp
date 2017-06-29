//FMExtendNode.cpp
#include "FMExtendNode.h"
#include "BWTAlgorithms.h"

FMExtendNode::FMExtendNode(const std::string* pQuery, FMExtendNode* parent) : m_pQuery(pQuery),m_pParent(parent)
{
}
// Destructor, recurisvely delete the children of the node
FMExtendNode::~FMExtendNode()
{
    // Delete children
    for(FMExtendNodePtrList::iterator iter = m_children.begin(); iter != m_children.end(); ++iter)
        delete *iter;
}


FMExtendNode* FMExtendNode::createChild(const std::string& lable,int direction)
{
		
	FMExtendNode* pAdded = new FMExtendNode(m_pQuery, this);
	m_children.push_back(pAdded);
	pAdded->extendToNewNode(lable,direction,left_label,right_label);
	return pAdded;
}
		
//Extend the label of this node by 1 bp , direction left:0 right:1
void FMExtendNode::extend(const std::string& extend_bp,int direction)
{
	if(direction==0)
		left_label.append(extend_bp);
	else if(direction==1)
		right_label.append(extend_bp);
}

void FMExtendNode::extendToNewNode(const std::string& extend_bp,int direction,std::string left,std::string right)
{
	if(direction==0)
	{
		left_label=left+extend_bp;
		right_label=right;
	}
	else if(direction==1)
	{
		left_label=left;
		right_label=right+extend_bp;
	}
	
}

//Set the initial kmer by two direction extension
void FMExtendNode::SetInitial()
{
	left_label="";
	right_label="";
}

std::pair<std::string,std::string> FMExtendNode::adjust_byInsertion(int direction,int seedSize,int InsertSize)
{
	//<preKmer,insertion_str>
	std::pair<std::string,std::string> out_strpair;
	std::string preKmer="";
	std::string insertion_str="";
	
	if(direction==0)
	{
		int insert_idx=(int)left_label.length()-seedSize;
		if(insert_idx>=0)
		{
		std::string pre_str=left_label.substr(0,insert_idx-InsertSize);
		std::string suffix_str=left_label.substr(insert_idx,seedSize);
		insertion_str=reverse(left_label.substr(insert_idx-InsertSize,InsertSize));
		left_label=pre_str+suffix_str;
		if((int)pre_str.length()>=seedSize)
			preKmer=reverse(pre_str.substr((int)pre_str.length()-seedSize,seedSize));
		else
			preKmer=reverse(pre_str);
		}
	}
	else
	{
		int insert_idx=(int)right_label.length()-seedSize;
		//printf("right_label-lth=%d\tinsert_idx=%d\tSeedSize=%d\tInsertSize=%d\n",(int)right_label.length(),insert_idx,seedSize,InsertSize);
		if(insert_idx>=0)
		{
		std::string pre_str=right_label.substr(0,insert_idx-InsertSize);
		std::string suffix_str=right_label.substr(insert_idx,seedSize);
		insertion_str=(right_label.substr(insert_idx-InsertSize,InsertSize));
		right_label=pre_str+suffix_str;
		if((int)pre_str.length()>=seedSize)
			preKmer=(pre_str.substr((int)pre_str.length()-seedSize,seedSize));
		else
			preKmer=(pre_str);
		}
	}
	
	out_strpair=std::make_pair(preKmer, insertion_str);
	return out_strpair;
}
void FMExtendNode::adjust_byDeletion(int direction,int seedSize,int deleSize)
{
	std::string dele_str="";
	for(int i=0;i<deleSize;i++)
		dele_str.append("-");
	if(direction==0)
	{
		int dele_idx=(int)left_label.length()-seedSize;
		if(dele_idx>=0)
		{
		std::string pre_str=left_label.substr(0,dele_idx);
		std::string suffix_str=left_label.substr(dele_idx,seedSize);
		left_label=pre_str+dele_str+suffix_str;
		}
	}
	else
	{
		int dele_idx=(int)right_label.length()-seedSize;
		if(dele_idx>=0)
		{
		std::string pre_str=right_label.substr(0,dele_idx);
		std::string suffix_str=right_label.substr(dele_idx,seedSize);
		right_label=pre_str+dele_str+suffix_str;
		}
	}
}
	
FMExOverlapNode* FMExOverlapNode::createChild(const std::string& label,int direction)
{
	FMExOverlapNode* PAdded = new FMExOverlapNode(m_pQuery,this);
	PAdded->extendToNewNode(label,direction,this->left_label,this->right_label);
	
	PAdded->lastSeedIdx=this->lastSeedIdx;//the newest match kmer index
	PAdded->nextSeedIdx=this->nextSeedIdx;//the next check kmer index on the QueryKmerArray;
	PAdded->numOfMatchSeeds=this->numOfMatchSeeds;
		
		//not necessary temporarily
		//int Left_start_idx;
		//int Right_end_idx;
		
	PAdded->initSeedIdx=this->initSeedIdx;//the initialSeed(kmer) idx on the Query
	PAdded->tempMismathCount=this->tempMismathCount;
	PAdded->insertion_count=this->insertion_count;
	PAdded->deletion_count=this->deletion_count;
	PAdded->mismatch_count=this->mismatch_count;
	
	PAdded->seed_rate=this->seed_rate;
	PAdded->isOverthreshold=this->isOverthreshold;
	//PAdded->mismatch_vector=this->mismatch_vector;
	//PAdded->deletion_vector=this->deletion_vector;
	PAdded->Insertion_vector=this->Insertion_vector;
	
	m_children.push_back(PAdded);
	
	return PAdded;
}






