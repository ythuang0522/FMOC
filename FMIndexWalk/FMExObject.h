#ifndef FMEXOBJECT_H
#define FMEXOBJECT_H

typedef std::pair<std::string,int> Insert_type;
typedef std::pair<std::string,int> Out_type;
typedef std::vector<Out_type> Out_test;
class Insertion
{
	public:
		int match_idx;//pre-kmer position on the Query
		int insert_size;
		std::string preKmer;
		std::string insertion_str;
		Insertion(int idx,int size,std::string pk,std::string insert_str)
		{
			match_idx=idx;
			insert_size=size;
			preKmer=pk;
			insertion_str=insert_str;
		}
		
};

class Solid_error
{
	public:
		int solid_left_idx;
		int solid_right_idx;
		std::string Left_error_1st_bp;
		std::string Right_error_1st_bp;
		Solid_error()
		{
		}
		Solid_error(int left_idx,int right_idx,std::string left_ebp,std::string right_ebp)
		{
			solid_left_idx=left_idx;
			solid_right_idx=right_idx;
			Left_error_1st_bp=left_ebp;
			Right_error_1st_bp=right_ebp;
		}
};
class Ksubstr
{
	public:
		std::string kmer;
		int countOfkmer;
		std::vector<Insert_type> In_str;
		
		//function
		Ksubstr(){}
		Ksubstr(std::string ksub,int count)
		{
			kmer=ksub;
			countOfkmer=count;
		}
		std::string getPreBp()
		{
			std::string out = kmer.substr(0,1);
			return out;
		}
		std::string getsufBp()
		{
			std::string out = kmer.substr((int)kmer.length()-1,1);
			return out;
		}
		std::string getPreStr()
		{
			std::string out = kmer.substr(0,(int)kmer.length()-1);
			return out;
		}
		std::string getsufStr()
		{
			std::string out = kmer.substr(1,(int)kmer.length()-1);
			return out;
		}
		
};

class OutInfo
{
	public:
		std::string outStr;
		int interval_size;
		int start_point;
		double Seed_rate;
		
		OutInfo()
		{
			outStr="";
			interval_size=0;
			start_point=-1;
			Seed_rate=0;
		}
		OutInfo(std::string inStr,int count,int idx,double rate)
		{
			outStr=inStr;
			interval_size=count;
			start_point=idx;
			Seed_rate=rate;
		}
};

typedef std::vector<Ksubstr> Ksub_vct;


#endif