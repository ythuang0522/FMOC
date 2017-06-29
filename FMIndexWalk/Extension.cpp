//Extension.cpp

#include "Extension.h"

typedef std::pair<int,int> start_end;

bool Extension::isLowComplexity (std::string seq, float threshold)
{
	size_t seqLen = seq.length();
	size_t countG =0 ;
	size_t countC =0 ;
	size_t countT =0 ;
	size_t countA =0 ;

	for (size_t i=0; i<seqLen; i++)
	{
		switch(seq[i]){
			case 'A': countA ++ ;break;
			case 'T': countT ++ ;break;
			case 'C': countC ++ ;break;
			case 'G': countG ++ ;break;
			default:  assert(false);
		}
	}
	if (  ((float) countA/seqLen >= threshold ) || ((float) countT/seqLen >= threshold)
			|| ((float) countC/seqLen >= threshold ) || ((float) countG/seqLen >= threshold) )
		return true;

	return false;
}

void Extension::getLRKmerInterval(std::string Query,int kmer_size,const BWT* pBWT, const BWT* pRBWT,
					std::vector<BWTInterval>& LTInterval,std::vector<BWTInterval>& RTInterval)
{
	/*
	for(int i = (int)(Query.length()-kmer_size); i >= 0 ; i -= 1)
	{
		std::string seedStr = Query.substr(i, kmer_size);
		LTInterval.push_back(BWTAlgorithms::findInterval( pBWT, seedStr) );		
	}
	*/
	for(int i = 0; i <= (int)(Query.length()-kmer_size) ; i += 1)
	{
		std::string seedStr = Query.substr(i, kmer_size);
		LTInterval.push_back(BWTAlgorithms::findInterval( pBWT, seedStr) );
		std::string seedStr2 = reverse( seedStr);
		RTInterval.push_back(BWTAlgorithms::findInterval( pRBWT, seedStr2) );	
	}
	/*
	for(int i = 0; i <= (int)(Query.length()-kmer_size) ; i += 1)
	{
		std::string seedStr = reverse( Query.substr(i, kmer_size) );
		RTInterval.push_back(BWTAlgorithms::findInterval( pRBWT, seedStr) );		
	}
	*/
}

void Extension::addQueryinOutInfo(std::string Query,std::vector<OutInfo>& Out_Info)
{
	OutInfo temp_outinfo(Query,1,0,1);
	Out_Info.push_back(temp_outinfo);
}
void Extension::addQueryinKsub(std::string Query,int kmer_size,std::vector<Ksub_vct>& correct_ksub)
{
	for(int k_idx=0;k_idx<=(int)Query.length()-kmer_size;k_idx++)
	{
		std::string ksub=Query.substr(k_idx,kmer_size);
		Ksubstr temp_ksub(ksub,1);
		Ksub_vct temp_ksub_vct;
		temp_ksub_vct.push_back(temp_ksub);
		correct_ksub.push_back(temp_ksub_vct);
	}
}

void Extension::addStrInKsub(std::vector<OutInfo>& Out_Info,int kmer_size,std::vector<Ksub_vct>& correct_ksub)
{
	for(int Size=0;Size<(int)Out_Info.size();Size++)
	{
		int end_idx=0;
		if( (int)correct_ksub.size()-Out_Info[Size].start_point < (int)Out_Info[Size].outStr.length()-kmer_size)
			end_idx=(int)correct_ksub.size()-Out_Info[Size].start_point;
		else
			end_idx=(int)Out_Info[Size].outStr.length()-kmer_size;
		
		for(int nowIdx=0;nowIdx<end_idx;nowIdx++)
		{
			/*
			if( nowIdx+Out_Info[Size].start_point >= (int)correct_ksub.size())
			{
				printf("Error with Out of range\t%d\t%d\tLength=\t%d\tStart=\t%d\n",nowIdx+Out_Info[Size].start_point,(int)correct_ksub.size(),
				(int)Out_Info[Size].outStr.length(),Out_Info[Size].start_point);
			}
			*/
			
			bool isExist=false;
			std::string temp_kmer=Out_Info[Size].outStr.substr(nowIdx,kmer_size);
			
			for(int t=0;t<(int)correct_ksub[nowIdx+Out_Info[Size].start_point].size();t++)
			{
				//printf("T_size=\t%d\n",t);
				if(correct_ksub[nowIdx+Out_Info[Size].start_point][t].kmer.compare(temp_kmer)==0)
				{
					isExist=true;
					correct_ksub[nowIdx+Out_Info[Size].start_point][t].countOfkmer+=Out_Info[Size].interval_size;
					break;
				}
			}
			
			if(!isExist)
			{
				Ksubstr temp_ksub(temp_kmer,Out_Info[Size].interval_size);
				//Ksub_vct temp_ksub_vct;
				correct_ksub[nowIdx+Out_Info[Size].start_point].push_back(temp_ksub);
			}
		}
	}
}

void Extension::addStrInKsub2(std::vector<OutInfo>& Out_Info,int kmer_size,std::vector<Out_test>& out_vct,std::string Query)
{
	printf("Count_of_out=\t%d\n",(int)Out_Info.size());
	for(int k_idx=0;k_idx<=(int)Query.length()-kmer_size;k_idx++)
	{
		Out_test temp_Out_test;
		for(int Size=0;Size<(int)Out_Info.size();Size++)
		{
			std::string temp_kmer;
			int m_idx=k_idx-Out_Info[Size].start_point;
			bool isExist=false;
			if( m_idx>=0 && (m_idx+kmer_size)<=(int)Out_Info[Size].outStr.length() )
			{
				temp_kmer=Out_Info[Size].outStr.substr(m_idx,kmer_size);
				for(int vct_size=0;vct_size<(int)temp_Out_test.size();vct_size++)
				{
					if(temp_Out_test[vct_size].first.compare(temp_kmer)==0)
					{
						isExist=true;
						break;
					}
				}
				if(!isExist)
					temp_Out_test.push_back(std::make_pair(temp_kmer,Out_Info[Size].interval_size));
			}
		}
		out_vct.push_back(temp_Out_test);
		
		//std::string ksub=Query.substr(k_idx,kmer_size);
		//Out_test temp_out;;
		//temp_out.push_back(std::make_pair(ksub,1));
		//out_vct.push_back(temp_out);
		//Ksubstr temp_ksub(ksub,1);
		//Ksub_vct temp_ksub_vct;
		//temp_ksub_vct.push_back(temp_ksub);
		//correct_ksub.push_back(temp_ksub_vct);
	}
	/*
	for(int Size=0;Size<(int)Out_Info.size();Size++)
	{
		int end_idx=0;
		if( ((int)Query.length()-kmer_size+1)-Out_Info[Size].start_point < (int)Out_Info[Size].outStr.length()-kmer_size)
			end_idx=((int)Query.length()-kmer_size+1)-Out_Info[Size].start_point;
		else
			end_idx=(int)Out_Info[Size].outStr.length()-kmer_size;
			
		for(int nowIdx=0;nowIdx<end_idx;nowIdx++)
		{
			bool isExist=false;
			std::string temp_kmer=Out_Info[Size].outStr.substr(nowIdx,kmer_size);
			
			for(int t=0;t<(int)out_vct[nowIdx+Out_Info[Size].start_point].size();t++)
			{
				//printf("T_size=\t%d\n",t);
				if(out_vct[nowIdx+Out_Info[Size].start_point][t].first.compare(temp_kmer)==0)
				{
					isExist=true;
					out_vct[nowIdx+Out_Info[Size].start_point][t].second+=Out_Info[Size].interval_size;
					break;
				}
			}
			
			if(!isExist)
			{
				//Ksubstr temp_ksub(temp_kmer,Out_Info[Size].interval_size);
				//Ksub_vct temp_ksub_vct;
				//correct_ksub[nowIdx+Out_Info[Size].start_point].push_back(temp_ksub);
				
				Out_test temp_out2;;
				temp_out2.push_back(std::make_pair(temp_kmer,Out_Info[Size].interval_size));
				out_vct.push_back(temp_out2);
			}
			
		}
	}
	*/
	
	
	/*
	for(int Size=0;Size<(int)Out_Info.size();Size++)
	{
		int end_idx=0;
		if( (int)correct_ksub.size()-Out_Info[Size].start_point < (int)Out_Info[Size].outStr.length()-kmer_size)
			end_idx=(int)correct_ksub.size()-Out_Info[Size].start_point;
		else
			end_idx=(int)Out_Info[Size].outStr.length()-kmer_size;
		
		for(int nowIdx=0;nowIdx<end_idx;nowIdx++)
		{
			
			bool isExist=false;
			std::string temp_kmer=Out_Info[Size].outStr.substr(nowIdx,kmer_size);
			
			for(int t=0;t<(int)correct_ksub[nowIdx+Out_Info[Size].start_point].size();t++)
			{
				//printf("T_size=\t%d\n",t);
				if(correct_ksub[nowIdx+Out_Info[Size].start_point][t].kmer.compare(temp_kmer)==0)
				{
					isExist=true;
					correct_ksub[nowIdx+Out_Info[Size].start_point][t].countOfkmer+=Out_Info[Size].interval_size;
					break;
				}
			}
			
			if(!isExist)
			{
				Ksubstr temp_ksub(temp_kmer,Out_Info[Size].interval_size);
				//Ksub_vct temp_ksub_vct;
				correct_ksub[nowIdx+Out_Info[Size].start_point].push_back(temp_ksub);
			}
			
		}
	}
	*/
}

void Extension::addStrInKsub3(std::vector<OutInfo>& Out_Info,int kmer_size,std::vector<Ksub_vct>& correct_ksub,std::string Query)
{
	//printf("Count_of_out=\t%d\n",(int)Out_Info.size());
	for(int k_idx=0;k_idx<=(int)Query.length()-kmer_size;k_idx++)
	{
		//Out_test temp_Out_test;
		Ksub_vct temp_ksub_vct;
		
		for(int Size=0;Size<(int)Out_Info.size();Size++)
		{
			std::string temp_kmer;
			int m_idx=k_idx-Out_Info[Size].start_point;
			bool isExist=false;
			if( m_idx>=0 && (m_idx+kmer_size)<=(int)Out_Info[Size].outStr.length() )
			{
				temp_kmer=Out_Info[Size].outStr.substr(m_idx,kmer_size);
				for(int vct_size=0;vct_size<(int)temp_ksub_vct.size();vct_size++)
				{
					if(temp_ksub_vct[vct_size].kmer.compare(temp_kmer)==0)
					{
						isExist=true;
						temp_ksub_vct[vct_size].countOfkmer+=Out_Info[Size].interval_size;
						break;
					}
				}
				if(!isExist)
				{
					Ksubstr in_ksub(temp_kmer,Out_Info[Size].interval_size);
					temp_ksub_vct.push_back(in_ksub);
				}
			}
		}
		correct_ksub.push_back(temp_ksub_vct);
		
		//std::string ksub=Query.substr(k_idx,kmer_size);
		//Out_test temp_out;;
		//temp_out.push_back(std::make_pair(ksub,1));
		//out_vct.push_back(temp_out);
		//Ksubstr temp_ksub(ksub,1);
		//Ksub_vct temp_ksub_vct;
		//temp_ksub_vct.push_back(temp_ksub);
		//correct_ksub.push_back(temp_ksub_vct);
	}
	
}

Solid_error Extension::getSolidRegion(std::string Query,int kmer_size,int thrshold,const BWT* pBWT)
{
	Solid_error out(-1,-1,"","");
	
	std::vector<BWTInterval> final_kmer_frq;
	std::vector<BWTInterval> final_kmer_frq_revc;
	std::vector<std::string> kmer_vct;
	//printf("K=%d\tQuery\t%s\n",kmer_size,Query.c_str());
	
	int start_idx=-1;
	int end_idx=-1;
	//std::pair<int,int> start_end;
	std::vector<start_end> s_e_vct;
	bool isSolid=false;
	bool isLargeThanOne=false;
	for(int i = 0; i <= (int)(Query.length())-kmer_size ; i++)
	{
		
			std::string seedStr = Query.substr(i, kmer_size);
			BWTInterval bip1=BWTAlgorithms::findInterval( pBWT, seedStr) ;
			BWTInterval bip1_revc=BWTAlgorithms::findInterval( pBWT, reverseComplement(seedStr)) ;
			final_kmer_frq.push_back(bip1);
			final_kmer_frq_revc.push_back(bip1_revc );
			kmer_vct.push_back(seedStr);
			
			if((int)bip1.size()>1 ||(int)bip1_revc.size()>1)
			{
				isLargeThanOne=true;
			}
			//if(!isSolid && bip1.isValid() && bip1_revc.isValid() &&( (int)bip1.size()>6 ||(int)bip1_revc.size()>6 ) )
			if(!isSolid && bip1.isValid() && bip1_revc.isValid() &&( (int)bip1.size()>thrshold ||(int)bip1_revc.size()>thrshold ) )
			{
				isSolid=true;
				start_idx=i;
				end_idx=i;
			}
			//else if( isSolid && bip1.isValid() && bip1_revc.isValid() &&( (int)bip1.size()>6 ||(int)bip1_revc.size()>6 ) )
			else if( isSolid && bip1.isValid() && bip1_revc.isValid() &&( (int)bip1.size()>thrshold ||(int)bip1_revc.size()>thrshold ) )
			{
				end_idx=i;
			}
			//else if( isSolid && (!bip1.isValid() || !bip1_revc.isValid()))
			//else if( isSolid && (!bip1.isValid() || !bip1_revc.isValid()||( (int)bip1.size()<=6 &&(int)bip1_revc.size()<=6 ) ))
			else if( isSolid && (!bip1.isValid() || !bip1_revc.isValid()||( (int)bip1.size()<=thrshold &&(int)bip1_revc.size()<=thrshold ) ))
			{
				s_e_vct.push_back(std::make_pair(start_idx,end_idx));
				isSolid=false;
			}
		
	}
/*
	for(int i=0;i<(int)final_kmer_frq.size();i++)
	{
		if(i>0)
		{
			// && (int)final_kmer_frq[i].size()>=5
			if((int)final_kmer_frq_revc[i-1].size()==0 && (int)final_kmer_frq_revc[i].size()>0)
			{
				printf("i-1=%d\tsize=%d\ti=%d\tsize=%d\n",i-1,(int)final_kmer_frq_revc[i-1].size(),i,(int)final_kmer_frq_revc[i].size());
				printf("Error at %d\terror_char=%s\tsuffix_kmer=%s\n",i-1,Query.substr(i-1,1).c_str(),kmer_vct[i].substr(0,kmer_size-1).c_str());
			}
			// && (int)final_kmer_frq[i-1].size()>=5
			else if( (int)final_kmer_frq_revc[i-1].size()>0 && (int)final_kmer_frq_revc[i].size()==0)
			{
				printf("i-1=%d\tsize=%d\ti=%d\tsize=%d\n",i-1,(int)final_kmer_frq_revc[i-1].size(),i,(int)final_kmer_frq_revc[i].size());
				printf("Error at %d\terror_char=%s\tpre_kmer=%s\n",i,kmer_vct[i].substr(kmer_size-1,1).c_str(),kmer_vct[i].substr(0,kmer_size-1).c_str());
			}
		}
	}
*/
	int max_solid_region_length=-1;
	for(int se_idx=0;se_idx<(int)s_e_vct.size();se_idx++)
	{
		int temp_length=s_e_vct[se_idx].second-s_e_vct[se_idx].first;
		if(temp_length>max_solid_region_length)
		{
			max_solid_region_length=temp_length;
			
			out.solid_left_idx=s_e_vct[se_idx].first;
			out.solid_right_idx=s_e_vct[se_idx].second;
			if(s_e_vct[se_idx].first>=1)
				out.Left_error_1st_bp=Query.substr(s_e_vct[se_idx].first-1,1);
			else
				out.Left_error_1st_bp="";
				
			if( s_e_vct[se_idx].second<(int)(Query.length())-kmer_size )
				out.Right_error_1st_bp=Query.substr(s_e_vct[se_idx].second+kmer_size,1);
			else
				out.Right_error_1st_bp="";
		}
		//printf("Start=\t%d\tEnd=\t%d\n",s_e_vct[se_idx].first,s_e_vct[se_idx].second);
	}
	/*
	printf("No.\t");
	for(int t=0;t<(int)final_kmer_frq.size();t++)
	{
		printf("%d\t",t);
	}
	printf("\nB\t");
	for(int t=0;t<(int)final_kmer_frq.size();t++)
	{
		printf("%d\t",(int)final_kmer_frq[t].size());
	}
	printf("\nB_r\t");
	for(int t=0;t<(int)final_kmer_frq.size();t++)
	{
		printf("%d\t",(int)final_kmer_frq_revc[t].size());
	}
	printf("\n");
	*/
	if( !isLargeThanOne )
		out.solid_left_idx=-3;
	
	
	return out;
}
Solid_error Extension::highError_getSolidRegion(std::string Query,int kmer_size,const BWT* pBWT)
{
	Solid_error out(-1,-1,"","");
	
	std::vector<BWTInterval> final_kmer_frq;
	std::vector<BWTInterval> final_kmer_frq_revc;
	//printf("K=%d\tQuery\t%s\n",kmer_size,Query.c_str());
	
	int start_idx=-1;
	int end_idx=-1;
	//std::pair<int,int> start_end;
	std::vector<start_end> s_e_vct;
	bool isSolid=false;
	
	for(int i = 0; i <= (int)(Query.length())-kmer_size ; i++)
	{
		
			std::string seedStr = Query.substr(i, kmer_size);
			BWTInterval bip1=BWTAlgorithms::findInterval( pBWT, seedStr) ;
			BWTInterval bip1_revc=BWTAlgorithms::findInterval( pBWT, reverseComplement(seedStr)) ;
			final_kmer_frq.push_back(bip1);
			final_kmer_frq_revc.push_back(bip1_revc );
			
			if(!isSolid &&( (int)bip1.size()>1 ||(int)bip1_revc.size()>1 ) )
			{
				isSolid=true;
				start_idx=i;
				end_idx=i;
			}
			else if( isSolid  &&( (int)bip1.size()>1 ||(int)bip1_revc.size()>1 ) )
			{
				end_idx=i;
			}
			//else if( isSolid && (!bip1.isValid() || !bip1_revc.isValid()))
			else if( isSolid && ( (int)bip1.size()<=1 &&(int)bip1_revc.size()<=1 ) )
			{
				s_e_vct.push_back(std::make_pair(start_idx,end_idx));
				isSolid=false;
			}
		
	}

	int max_solid_region_length=-1;
	for(int se_idx=0;se_idx<(int)s_e_vct.size();se_idx++)
	{
		int temp_length=s_e_vct[se_idx].second-s_e_vct[se_idx].first;
		if(temp_length>max_solid_region_length)
		{
			max_solid_region_length=temp_length;
			
			out.solid_left_idx=s_e_vct[se_idx].first;
			out.solid_right_idx=s_e_vct[se_idx].second;
			if(s_e_vct[se_idx].first>=1)
				out.Left_error_1st_bp=Query.substr(s_e_vct[se_idx].first-1,1);
			else
				out.Left_error_1st_bp="";
				
			if( s_e_vct[se_idx].second<(int)(Query.length())-kmer_size )
				out.Right_error_1st_bp=Query.substr(s_e_vct[se_idx].second+kmer_size,1);
			else
				out.Right_error_1st_bp="";
		}
		//printf("Start=\t%d\tEnd=\t%d\n",s_e_vct[se_idx].first,s_e_vct[se_idx].second);
	}
	
	
	return out;
}
std::string Extension::getSolidRegion_v2(std::string Query,int Seed_size,const BWT* pBWT)
{
	std::vector<BWTInterval> final_kmer_frq;
	std::vector<BWTInterval> final_kmer_frq_revc;
	std::string consensus="";
	int start_idx=-1;
	int end_idx=-1;
	bool have_highFQ=false;
	std::vector<int> highFQ;
	
	for(int i = 0; i <= (int)(Query.length()-Seed_size) ; i++)
	{
		std::string seedStr = Query.substr(i, Seed_size);
		BWTInterval bip1=BWTAlgorithms::findInterval( pBWT, seedStr) ;
		BWTInterval bip1_revc=BWTAlgorithms::findInterval( pBWT, reverseComplement(seedStr)) ;
		
		final_kmer_frq.push_back(bip1);
		final_kmer_frq_revc.push_back(bip1_revc );
		
		if( ((int)bip1.size()>=5 && (int)bip1_revc.size() >=1) || ((int)bip1.size()>=1 && (int)bip1_revc.size() >=5) )
		{
			if( ((int)bip1.size()>=52&& (int)final_kmer_frq[i-1].size()<52)||((int)bip1_revc.size() >=52 && (int)final_kmer_frq_revc[i-1].size()<52) )
			{
				highFQ.push_back(i);
				have_highFQ=true;
			}
			
			if(start_idx>=0)
			{
				end_idx=i;
			}
			else{
				
				start_idx=i;
				end_idx=i;
			}
		}
		else{
			if(start_idx>=0)
			{
				break;
			}
			
		}
	}
	if(start_idx>=0 && end_idx>=0)
	{
		if(have_highFQ)
		{
			if( (int)highFQ.size()>=2)
			{
				printf("Query=\t%s\nSolid_left_idx=%d\tSolid_right_idx=%d\n",Query.c_str(),start_idx,end_idx);
				consensus=Query.substr(start_idx,(int)Seed_size+(end_idx-start_idx));
			}
			else if((int)highFQ.size()==1)
			{
				printf("Query=\t%s\nSolid_left_idx=%d\tSolid_right_idx=%d\n",Query.c_str(),highFQ[0],end_idx);
				consensus=Query.substr(highFQ[0],(int)Seed_size+(end_idx-highFQ[0]));
			}
		}
		else{
		printf("Query=\t%s\nSolid_left_idx=%d\tSolid_right_idx=%d\n",Query.c_str(),start_idx,end_idx);
		consensus=Query.substr(start_idx,(int)Seed_size+(end_idx-start_idx));
		}
	}
	else
	{
		consensus=Query;
	}
	printf("Correct=\t%s\n",consensus.c_str());
	printf("No.\t");
	for(int t=0;t<(int)final_kmer_frq.size();t++)
	{
		printf("%d\t",t);
	}
	printf("\nB\t");
	for(int t=0;t<(int)final_kmer_frq.size();t++)
	{
		printf("%d\t",(int)final_kmer_frq[t].size());
	}
	printf("\nB_r\t");
	for(int t=0;t<(int)final_kmer_frq.size();t++)
	{
		printf("%d\t",(int)final_kmer_frq_revc[t].size());
	}
	printf("\n");
	return consensus;
}

bool Extension::ExtensionRead(std::string Query,int kmer_size,int check_kmer_size,int solid_idx,const BWT* pBWT, const BWT* pRBWT,
					std::vector<OutInfo>& Out_Info,std::vector<BWTInterval>& LTInterval,std::vector<BWTInterval>& RTInterval)
{
	int temp_idx=solid_idx;
	solid_idx+=1;
	BWTIntervalPair bip;
	std::string initSeed;
	do
	{
		solid_idx-=1;
		if(solid_idx>=0)
		{
			initSeed = Query.substr(solid_idx,kmer_size);
			//printf("The Seed kmer = %s\n",initSeed.c_str());
			bip = BWTAlgorithms::findIntervalPair( pBWT, pRBWT, initSeed);
		}
	}while((int)bip.interval[0].size() >60 && solid_idx>=0);
	
	if(solid_idx<0)
	{
		solid_idx=temp_idx;
		initSeed = Query.substr(solid_idx,kmer_size);
		bip = BWTAlgorithms::findIntervalPair( pBWT, pRBWT, initSeed);
	}
	//if( bip.isValid() && (bip.interval[0].size() >=10 || bip.interval[1].size() >=10))
	//{
		//if(!isLowComplexity(initSeed))
		//{
			//printf("Kmer=\t%s\tSize1=\t%d\tSize2=\t%d\n",initSeed.c_str(),(int)bip.interval[0].size(),(int)bip.interval[1].size());
			FMExtendTree OverlapTree(Query,initSeed,check_kmer_size,bip,solid_idx,(int)Query.length()/2,0.9,pBWT,pRBWT,Out_Info,LTInterval,RTInterval);
			
			//printf("Left-extension start.\n");
			while(OverlapTree.getCurrentLength_L() < Query.length()+20)
			{
				if(OverlapTree.isEmpty()) break;
			
				int flag = OverlapTree.extendOneBase(0);
				
				if(flag == -3)
				{
					//printf("The_leaves_are_too_much\tThe kmer_interval=\t%d\n",(int)bip.interval[0].size() );
					break;
				}
				if(flag == -1)
				{
					//printf("GOOD_JOB\n");
					break;
				}
				if(flag == -2)
				{
					//printf("exceed search depth\n");
					break;
				}
			}
			//printf("Left-extension end.\n");
			
			OverlapTree.Left2Right();
			//printf("Right-extension start.\n");
			
			while(OverlapTree.getCurrentLength_R() < Query.length()+20)
			{
				if(OverlapTree.isEmpty()) break;
				
				int flag = OverlapTree.extendOneBase(1);
				if(flag == -3)
				{
					//printf("The_leaves_are_too_much\tThe kmer_interval=\t%d\n",(int)bip.interval[0].size() );
					break;
				}
				if(flag == -1)
				{
					//printf("high error\n");
					break;
				}
				if(flag == -2)
				{
					//printf("exceed search depth\n");
					break;
				}
			}
			
			//printf("Right-extension end.\n");
			
	//}
	
	std::string RC_initSeed = reverseComplement(initSeed);
	BWTIntervalPair RC_bip = BWTAlgorithms::findIntervalPair( pBWT, pRBWT, RC_initSeed);
	
	FMExtendTreeRC OverlapTreeRC(Query,RC_initSeed,check_kmer_size,RC_bip,solid_idx,(int)Query.length()/2,0.9,pBWT,pRBWT,Out_Info,LTInterval,RTInterval);
			
			//printf("Left-extension start.\n");
			while(OverlapTreeRC.getCurrentLength_L() < Query.length()+20)
			{
				if(OverlapTreeRC.isEmpty()) break;
			
				int flag = OverlapTreeRC.extendOneBase(0);
				
				if(flag == -3)
				{
					//printf("The_leaves_are_too_much\tThe kmer_interval=\t%d\n",(int)bip.interval[0].size() );
					break;
				}
				if(flag == -1)
				{
					//printf("GOOD_JOB\n");
					break;
				}
				if(flag == -2)
				{
					//printf("exceed search depth\n");
					break;
				}
			}
			//printf("Left-extension end.\n");
			
			OverlapTreeRC.Left2Right();
			//printf("Right-extension start.\n");
			
			while(OverlapTreeRC.getCurrentLength_R() < Query.length()+20)
			{
				if(OverlapTreeRC.isEmpty()) break;
				
				int flag = OverlapTreeRC.extendOneBase(1);
				if(flag == -3)
				{
					//printf("The_leaves_are_too_much\tThe kmer_interval=\t%d\n",(int)bip.interval[0].size() );
					break;
				}
				if(flag == -1)
				{
					//printf("high error\n");
					break;
				}
				if(flag == -2)
				{
					//printf("exceed search depth\n");
					break;
				}
			}
	return true;
}

std::string Extension::correct_right(Solid_error& solid_info,std::vector<Ksub_vct>& correct_ksub,int k_diff,std::string last_kmer)
{
	std::string right_out="";
	std::string curr_kmer=last_kmer;
	int start_idx=solid_info.solid_right_idx+k_diff+1;
	bool isFirst=false;
	
	for(int i=start_idx;i<(int)correct_ksub.size();i++)
	{
		int max_idx=-1;
		int max_count=-1;
		if(isFirst)
		{
			for(int kmer_idx=0;kmer_idx<(int)correct_ksub[i].size();kmer_idx++)
			{
				if( correct_ksub[i][kmer_idx].getPreStr().compare(curr_kmer)==0 )
				{
					if( correct_ksub[i][kmer_idx].getsufBp().compare(solid_info.Right_error_1st_bp)!=0 )
					{
						if(max_count<correct_ksub[i][kmer_idx].countOfkmer)
						{
							max_idx=kmer_idx;
							max_count=correct_ksub[i][kmer_idx].countOfkmer;
						}
					}
				}
			}
			if(max_idx>=0)
			{
				if(correct_ksub[i][max_idx].getsufBp().compare("-")!=0)
				{
					right_out.append(correct_ksub[i][max_idx].getsufBp());
				}
				curr_kmer=correct_ksub[i][max_idx].getsufStr();
			}
			else
			{
				break;
			}
			isFirst=false;
		}
		else
		{
			for(int kmer_idx=0;kmer_idx<(int)correct_ksub[i].size();kmer_idx++)
			{
				if( correct_ksub[i][kmer_idx].getPreStr().compare(curr_kmer)==0 )
				{
					if(max_count<correct_ksub[i][kmer_idx].countOfkmer)
					{
						max_idx=kmer_idx;
						max_count=correct_ksub[i][kmer_idx].countOfkmer;
					}
				}
			}
			if(max_idx>=0)
			{
				if(correct_ksub[i][max_idx].getsufBp().compare("-")!=0)
				{
					right_out.append(correct_ksub[i][max_idx].getsufBp());
				}
				curr_kmer=correct_ksub[i][max_idx].getsufStr();
			}
			else
			{
				break;
			}
		}
	}
	
	return right_out;
}

std::string Extension::correct_left(Solid_error& solid_info,std::vector<Ksub_vct>& correct_ksub,std::string last_kmer)
{
	std::string left_out="";
	std::string curr_kmer=last_kmer;
	int start_idx=solid_info.solid_left_idx-1;
	bool isFirst=false;
	
	for(int i=start_idx;i>=0;i--)
	{
		int max_idx=-1;
		int max_count=-1;
		if(isFirst)
		{
			for(int kmer_idx=0;kmer_idx<(int)correct_ksub[i].size();kmer_idx++)
			{
				if( correct_ksub[i][kmer_idx].getsufStr().compare(curr_kmer)==0 )
				{
					if( correct_ksub[i][kmer_idx].getPreBp().compare(solid_info.Left_error_1st_bp)!=0 )
					{
						if(max_count<correct_ksub[i][kmer_idx].countOfkmer)
						{
							max_idx=kmer_idx;
							max_count=correct_ksub[i][kmer_idx].countOfkmer;
						}
					}
				}
			}
			if(max_idx>=0)
			{
				if(correct_ksub[i][max_idx].getPreBp().compare("-")!=0)
				{
					left_out.append(correct_ksub[i][max_idx].getPreBp());
				}
				curr_kmer=correct_ksub[i][max_idx].getPreStr();
			}
			else
			{
				break;
			}
			isFirst=false;
		}
		else
		{
			for(int kmer_idx=0;kmer_idx<(int)correct_ksub[i].size();kmer_idx++)
			{
				if( correct_ksub[i][kmer_idx].getsufStr().compare(curr_kmer)==0 )
				{
					if(max_count<correct_ksub[i][kmer_idx].countOfkmer)
					{
						max_idx=kmer_idx;
						max_count=correct_ksub[i][kmer_idx].countOfkmer;
					}
				}
			}
			if(max_idx>=0)
			{
				if(correct_ksub[i][max_idx].getPreBp().compare("-")!=0)
				{
					left_out=correct_ksub[i][max_idx].getPreBp()+left_out;
					//left_out.append(correct_ksub[i][max_idx].getPreBp());
				}
				curr_kmer=correct_ksub[i][max_idx].getPreStr();
			}
			else
			{
				break;
			}
		}
	}
	return left_out;
	//return reverse(left_out);
}

std::string Extension::Newcorrect_right(std::string Query,Solid_error& solid_info,std::vector<Ksub_vct>& correct_ksub,int k_diff,std::string last_kmer,std::string FirstKmer,int Solid_kmer_size,const BWT* pBWT)
{
	std::string right_out=FirstKmer;
	std::string true_right_out="";
	std::string curr_kmer=last_kmer;
	int start_idx=solid_info.solid_right_idx+k_diff+1;
	std::vector<int> curr_MaxIdx;
	curr_MaxIdx.clear();
	std::string no_crStr="";
	
	for(int i=start_idx;i<(int)correct_ksub.size();i++)
	{
		int max_idx=-1;
		int max_count=-1;
		
		
		for(int kmer_idx=0;kmer_idx<(int)correct_ksub[i].size();kmer_idx++)
		{
			if( correct_ksub[i][kmer_idx].getPreStr().compare(curr_kmer)==0 )
			{
				if(max_count<correct_ksub[i][kmer_idx].countOfkmer)
				{
					max_idx=kmer_idx;
					max_count=correct_ksub[i][kmer_idx].countOfkmer;
					curr_MaxIdx.clear();
					curr_MaxIdx.push_back(kmer_idx);
				}
				else if(max_count==correct_ksub[i][kmer_idx].countOfkmer)
				{
					curr_MaxIdx.push_back(kmer_idx);
				}
			}
		}
		if(max_idx>=0)
		{
			if((int)curr_MaxIdx.size() > 1 )
			{
				max_idx=-1;
				max_count=-1;
				
				for(int t=0;t<(int)curr_MaxIdx.size();t++)
				{
					std::string check_kmer="";
					int temp_max_count=correct_ksub[i][curr_MaxIdx[t]].countOfkmer;;
					if(correct_ksub[i][curr_MaxIdx[t]].getsufBp().compare("-")!=0)
					{
						check_kmer=right_out.substr((int)right_out.length()-(Solid_kmer_size-1),Solid_kmer_size-1)+correct_ksub[i][curr_MaxIdx[t]].getsufBp();
						
						BWTInterval bip1_revc=BWTAlgorithms::findInterval(pBWT, reverseComplement(check_kmer));
						if( bip1_revc.isValid())
						{
							//printf("At %d:\t%s\t%d,%d\n",i,check_kmer.c_str(),temp_max_count,(int)bip1_revc.size());	
							temp_max_count+=(int)bip1_revc.size();
							
						}	
					}
					
					if( temp_max_count > max_count)
					{
						max_idx=curr_MaxIdx[t];
						max_count=temp_max_count;
					}
					
				}
			}
			
			if(correct_ksub[i][max_idx].getsufBp().compare("-")!=0)
			{
				right_out.append(correct_ksub[i][max_idx].getsufBp());
				true_right_out.append(correct_ksub[i][max_idx].getsufBp());
			}
			curr_kmer=correct_ksub[i][max_idx].getsufStr();
		}
		else
		{
			int Q_length=(int)Query.length();
			int no_crIdx=i+(-k_diff+Solid_kmer_size-1);
			no_crStr=Query.substr(no_crIdx,Q_length-no_crIdx);
			break;
		}
		
	}
	
	//return true_right_out+no_crStr;
	return true_right_out;
}

std::string Extension::Newcorrect_left(std::string Query,Solid_error& solid_info,std::vector<Ksub_vct>& correct_ksub,std::string last_kmer,std::string FirstKmer,int Solid_kmer_size,const BWT* pBWT)
{
	std::string left_out=FirstKmer;
	std::string true_left_out="";
	std::string curr_kmer=last_kmer;
	int start_idx=solid_info.solid_left_idx-1;
	std::vector<int> curr_MaxIdx;
	curr_MaxIdx.clear();
	std::string no_crStr="";
	
	for(int i=start_idx;i>=0;i--)
	{
		int max_idx=-1;
		int max_count=-1;
		
		
		for(int kmer_idx=0;kmer_idx<(int)correct_ksub[i].size();kmer_idx++)
		{
			if( correct_ksub[i][kmer_idx].getsufStr().compare(curr_kmer)==0 )
			{
				if(max_count<correct_ksub[i][kmer_idx].countOfkmer)
				{
					max_idx=kmer_idx;
					max_count=correct_ksub[i][kmer_idx].countOfkmer;
					curr_MaxIdx.clear();
					curr_MaxIdx.push_back(kmer_idx);
				}
				else if(max_count==correct_ksub[i][kmer_idx].countOfkmer)
				{
					curr_MaxIdx.push_back(kmer_idx);
				}
			}
		}
		
		
		if(max_idx>=0)
		{
			if((int)curr_MaxIdx.size() > 1 )
			{
				max_idx=-1;
				max_count=-1;
				
				for(int t=0;t<(int)curr_MaxIdx.size();t++)
				{
					std::string check_kmer="";
					int temp_max_count=correct_ksub[i][curr_MaxIdx[t]].countOfkmer;;
					if(correct_ksub[i][curr_MaxIdx[t]].getPreBp().compare("-")!=0)
					{
						check_kmer=correct_ksub[i][curr_MaxIdx[t]].getPreBp()+left_out.substr(0,Solid_kmer_size-1);
						
						BWTInterval bip1=BWTAlgorithms::findInterval( pBWT, check_kmer);
						BWTInterval bip1_revc=BWTAlgorithms::findInterval( pBWT, reverseComplement(check_kmer));
						if( bip1.isValid() && bip1_revc.isValid())
						{
							//printf("At %d:\t%s\t%d,%d\n",i,check_kmer.c_str(),temp_max_count,(int)bip1_revc.size());	
							temp_max_count+=(int)bip1_revc.size();
							
						}	
					}
					
					if( temp_max_count > max_count)
					{
						max_idx=curr_MaxIdx[t];
						max_count=temp_max_count;
					}
					
				}
			}
		
		
		
		
			if(correct_ksub[i][max_idx].getPreBp().compare("-")!=0)
			{
				left_out=correct_ksub[i][max_idx].getPreBp()+left_out;
				true_left_out.append(correct_ksub[i][max_idx].getPreBp());
				//left_out.append(correct_ksub[i][max_idx].getPreBp());
			}
			curr_kmer=correct_ksub[i][max_idx].getPreStr();
		}
		else
		{
			no_crStr=Query.substr(0,i);
			break;
		}
		
	}
	//return left_out;
	//return no_crStr+reverse(true_left_out);
	return reverse(true_left_out);
}

std::string Extension::TrimReads(std::string consensus,int kmer_size,const BWT* pBWT)
{
	//printf("CST=\t%s\n",consensus.c_str());
	std::string out="";
	
	std::vector<BWTInterval> final_kmer_frq;
	std::vector<BWTInterval> final_kmer_frq_revc;
	
	int start_idx=-1;
	int end_idx=-1;
	//std::pair<int,int> start_end;
	std::vector<start_end> s_e_vct;
	bool isSolid=false;
	
	for(int i = 0; i <= (int)(consensus.length())-kmer_size ; i++)
	{
			std::string seedStr = consensus.substr(i, kmer_size);
			BWTInterval bip1=BWTAlgorithms::findInterval( pBWT, seedStr) ;
			BWTInterval bip1_revc=BWTAlgorithms::findInterval( pBWT, reverseComplement(seedStr)) ;
			final_kmer_frq.push_back(bip1);
			final_kmer_frq_revc.push_back(bip1_revc );
			if(!isSolid && bip1.isValid() && bip1_revc.isValid() &&( (int)bip1.size()>6 ||(int)bip1_revc.size()>6 ) )
			{
				isSolid=true;
				start_idx=i;
				end_idx=i;
			}
			else if( isSolid && bip1.isValid() && bip1_revc.isValid() &&( (int)bip1.size()>6 ||(int)bip1_revc.size()>6 ) )
			{
				end_idx=i;
			}
			else if( isSolid && (!bip1.isValid() || !bip1_revc.isValid()||( (int)bip1.size()<=6 &&(int)bip1_revc.size()<=6 ) ))
			{
				s_e_vct.push_back(std::make_pair(start_idx,end_idx));
				isSolid=false;
			}
		
	}

	//int max_solid_region_length=-1;
	int left_idx=-1;
	int right_idx=-1;
	for(int se_idx=0;se_idx<(int)s_e_vct.size();se_idx++)
	{
		if(se_idx==0)
		{
			left_idx=s_e_vct[se_idx].first;
			right_idx=s_e_vct[se_idx].second;
		}
		else
		{
			if(right_idx<s_e_vct[se_idx].second)
			{
				right_idx=s_e_vct[se_idx].second;
			}
		}
		/*
		int temp_length=s_e_vct[se_idx].second-s_e_vct[se_idx].first;
		if(temp_length>max_solid_region_length)
		{
			max_solid_region_length=temp_length;
			
			left_idx=s_e_vct[se_idx].first;
			right_idx=s_e_vct[se_idx].second;
			
		}
		*/
		//printf("Start=\t%d\tEnd=\t%d\n",s_e_vct[se_idx].first,s_e_vct[se_idx].second);
	}
	if(left_idx>=0 && right_idx>=0)
		out=consensus.substr(left_idx,right_idx-left_idx+kmer_size);
	else
		out=consensus;
	
	printf("No.\t");
	for(int t=0;t<(int)final_kmer_frq.size();t++)
	{
		printf("%d\t",t);
	}
	printf("\nB\t");
	for(int t=0;t<(int)final_kmer_frq.size();t++)
	{
		printf("%d\t",(int)final_kmer_frq[t].size());
	}
	printf("\nB_r\t");
	for(int t=0;t<(int)final_kmer_frq.size();t++)
	{
		printf("%d\t",(int)final_kmer_frq_revc[t].size());
	}
	printf("\n");
	
	return out;
}

std::string Extension::NewTrimReads(std::string consensus,int kmer_size,const BWT* pBWT)
{
	std::string out="";
	
	std::vector<BWTInterval> final_kmer_frq;
	std::vector<BWTInterval> final_kmer_frq_revc;
	
	int start_idx=-1;
	int end_idx=-1;
	
	//Left->Right
	for(int i = 0; i <= (int)(consensus.length())-kmer_size ; i++)
	{
		std::string seedStr = consensus.substr(i, kmer_size);
		BWTInterval bip1=BWTAlgorithms::findInterval( pBWT, seedStr) ;
		BWTInterval bip1_revc=BWTAlgorithms::findInterval( pBWT, reverseComplement(seedStr)) ;
		final_kmer_frq.push_back(bip1);
		final_kmer_frq_revc.push_back(bip1_revc );
		if( bip1.isValid() && bip1_revc.isValid() &&( (int)bip1.size()>6 ||(int)bip1_revc.size()>6 ) )
		{
			start_idx=i;
			break;
		}
	}
	//Right->Left
	for(int i = (int)(consensus.length())-kmer_size; i >= 0 ; i--)
	{
		std::string seedStr = consensus.substr(i, kmer_size);
		BWTInterval bip1=BWTAlgorithms::findInterval( pBWT, seedStr) ;
		BWTInterval bip1_revc=BWTAlgorithms::findInterval( pBWT, reverseComplement(seedStr)) ;
		final_kmer_frq.push_back(bip1);
		final_kmer_frq_revc.push_back(bip1_revc );
		if( bip1.isValid() && bip1_revc.isValid() &&( (int)bip1.size()>6 ||(int)bip1_revc.size()>6 ) )
		{
			end_idx=i;
			break;
		}
	}
	
	
	if(start_idx>=0 && end_idx>=0 && (end_idx>=start_idx))
		out=consensus.substr(start_idx,end_idx-start_idx+kmer_size);
	else
		out=consensus;
	
	
	return out;
}

void Extension::printfKFQ(std::string Query,int kmer_size,const BWT* pBWT)
{
	std::vector<BWTInterval> final_kmer_frq;
	std::vector<BWTInterval> final_kmer_frq_revc;
	//std::vector<BWTInterval> final_kmer_frq_complement;
	
	for(int i = 0; i <= (int)(Query.length())-kmer_size ; i++)
	{
		
			std::string seedStr = Query.substr(i, kmer_size);
			BWTInterval bip1=BWTAlgorithms::findInterval( pBWT, seedStr) ;
			BWTInterval bip1_revc=BWTAlgorithms::findInterval( pBWT, reverseComplement(seedStr)) ;
			//BWTInterval bip1_complement=BWTAlgorithms::findInterval( pBWT, complement(seedStr)) ;
			final_kmer_frq.push_back(bip1);
			final_kmer_frq_revc.push_back(bip1_revc );
			//final_kmer_frq_complement.push_back(bip1_complement );
		
	}
	
	printf("No.\t");
	for(int t=0;t<(int)final_kmer_frq.size();t++)
	{
		printf("%d\t",t);
	}
	printf("\nB\t");
	for(int t=0;t<(int)final_kmer_frq.size();t++)
	{
		printf("%d\t",(int)final_kmer_frq[t].size());
	}
	printf("\nB_r\t");
	for(int t=0;t<(int)final_kmer_frq.size();t++)
	{
		printf("%d\t",(int)final_kmer_frq_revc[t].size());
	}
	/*
	printf("\nB_c\t");
	for(int t=0;t<(int)final_kmer_frq.size();t++)
	{
		printf("%d\t",(int)final_kmer_frq_complement[t].size());
	}
	*/
	printf("\n");
	

}