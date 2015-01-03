/*
 * GeneInteractionBinary.cpp
 *
 *  Created on: Jan 1, 2015
 *      Author: ionadmin
 */

#include "GeneInteractionBinary.h"
#include <algorithm>


namespace SHEsis {





GeneInteractionBinary::GeneInteractionBinary(boost::shared_ptr<SHEsisData> data):GeneInteraction(data) {
	// TODO Auto-generated constructor stub

}

GeneInteractionBinary::~GeneInteractionBinary() {
	// TODO Auto-generated destructor stub
}

void GeneInteractionBinary::print(){
	for(int i=0;i<this->res.size();i++){
		std::cout<<res[i].snpset<<"\tcontrol:"<<res[i].ctrlEntropy<<"\tcase:"<<res[i].caseEntropy<<"\tdiff:"<<res[i].diff<<"\n";
	}
}
void printvector2D(std::vector<std::vector<int> > snp){
	for(int i=0;i<snp.size();i++){
		for(int j=0;j<snp[i].size();j++){
			std::cout<<snp[i][j];
		}
		std::cout<<"\n";
	}
}
void printvector1D(std::vector<int> snp){
	for(int i=0;i<snp.size();i++){
			std::cout<<snp[i]<<",";
	}
	std::cout<<"\n";
}

void GeneInteractionBinary::CalGeneInteraction(){
	BOOST_ASSERT(this->lowbound<=this->upperbound);
	for(int i=this->lowbound;i<=this->upperbound;i++){
		std::vector<std::vector<int>  > Snp;
		GenerateSNPCombination(i,Snp);
		for(int j=0;j<Snp.size();j++){
			this->res.push_back(this->GetOneSNPCombinationInformationGain(Snp[j]));
		}
	}
}

gxgBinaryRes GeneInteractionBinary::GetOneSNPCombinationInformationGain(std::vector<int>& Snp){
	gxgBinaryRes ret;
	std::string snpset="";
	for(int i=0;i<Snp.size()-1;i++){
		snpset+=this->data->vLocusName[Snp[i]]+",";
	}
	snpset+=this->data->vLocusName[Snp[Snp.size()-1]];
	ret.snpset=snpset;
//	std::cout<<"calculating "<<snpset<<"\n";
	double caseEntropy=0,controlEntropy=0;
	std::vector<std::vector<std::string> > cp;
	this->GenerateGenotypeCombination(Snp,cp);
	this->GetInformationGain(Snp,cp,ret.caseEntropy,ret.ctrlEntropy);
	ret.diff=ret.caseEntropy-ret.ctrlEntropy;
	return ret;
}

void GeneInteractionBinary::GenerateSNPCombination(int snpnum,std::vector<std::vector<int> >& ret){
	ret.clear();
	std::vector<int> r;
	std::vector<int> c;

	BOOST_ASSERT(snpnum<=this->data->getSnpNum());
	for(int i=0; i <snpnum;i++)
		r.push_back(i);
//	std::cout<<"snpnum:"<<this->data->getSnpNum()<<"\n";
	for(int i=0;i<this->data->getSnpNum();i++)
		c.push_back(i);

	  do {
	    std::vector<int> tmp = r;
	    ret.push_back(tmp);
	  } while (next_combination(c.begin(), c.end(), r.begin(), r.end()));
}

void GeneInteractionBinary::GetInformationGain(std::vector<int>& Snp,std::vector<std::vector<std::string> >& cp,double& caseGain,double& ctrlGain)
{
	caseGain=0;
	ctrlGain=0;
	std::vector<int> validCase;
	std::vector<int> validCtrl;
	double caseSum,ctrlSum;
	this->getNonmissingSample(Snp,validCase,validCtrl);
//	std::cout<<"valid case:";
//	printvector1D(validCase);
//	std::cout<<"valid control:";
//	printvector1D(validCtrl);
	this->getSingleEntropySum(Snp,validCase,validCtrl,caseSum,ctrlSum);
//	std::cout<<"case Entropy sum="<<caseSum<<"\n";
//	std::cout<<"ctrl Entropy sum="<<ctrlSum<<"\n";
	//calculate case information gain
	double caseMutual=0;
	double ctrlMutual=0;
	for(int cpIdx=0;cpIdx<cp.size();cpIdx++){
		int count=0;
		for(int sample=0;sample<validCase.size();sample++){
			bool equal=true;
			for(int i=0;i<cp[cpIdx].size();i++){
				if(EQUAL != this->genotypeEqual(cp[cpIdx][i],validCase[sample],Snp[i])){
					equal = false;
					break;
				}
			}
			if(equal)
				count++;
		}
		double rate=(double)count/(double)validCase.size();
		if(rate!=0)
			caseMutual+=rate*log(rate);
	}

	//calculate case information gain
	for(int cpIdx=0;cpIdx<cp.size();cpIdx++){
		int count=0;
		for(int sample=0;sample<validCtrl.size();sample++){
			bool equal=true;
			for(int i=0;i<cp[cpIdx].size();i++){
				if(EQUAL != this->genotypeEqual(cp[cpIdx][i],validCtrl[sample],Snp[i])){
					equal = false;
					break;
				}
			}
			if(equal)
				count++;
		}
		double rate=(double)count/(double)validCtrl.size();
		if(rate!=0)
			ctrlMutual+=rate*log(rate);
	}
	caseMutual*=(-1);
	ctrlMutual*=(-1);
//	std::cout<<"case Mutual entropy="<<caseMutual<<"\n";
//	std::cout<<"ctrl Mutual entropy="<<ctrlMutual<<"\n";
	caseGain=caseSum-caseMutual;
	ctrlGain=ctrlSum-ctrlMutual;
}

void GeneInteractionBinary::getNonmissingSample(std::vector<int>& Snp,std::vector<int>& validCase,std::vector<int>& validCtrl){
	validCase.clear();
	validCtrl.clear();
	for(int sample=0;sample<this->data->getSampleNum();sample++){
		bool valid=true;
		for(int i=0;i<Snp.size();i++){
			if(this->data->mGenotype[sample][Snp[i]][0]==0){
				valid=false;
				break;
			}
		}
		if(valid){
			if(this->data->vLabel[sample] == CASE)
				validCase.push_back(sample);
			else if (this->data->vLabel[sample] == CONTROL)
				validCtrl.push_back(sample);
		}
	}
}

void GeneInteractionBinary::getSingleEntropySum(std::vector<int>& Snp,std::vector<int>& validCase,std::vector<int>& validCtrl,double& CaseSum, double& CtrlSum)
{
	CaseSum=0;
	CtrlSum=0;
	boost::unordered_map<std::string,int> caseCount;
	boost::unordered_map<std::string,int> ctrlCount;
	for(int snp=0;snp<Snp.size();snp++){
		caseCount.clear();
		ctrlCount.clear();
		for(int _case=0;_case<validCase.size();_case++){
			if(caseCount.find(this->mGenotypeStr[validCase[_case]][Snp[snp]]) != caseCount.end()){
				caseCount[this->mGenotypeStr[validCase[_case]][Snp[snp]]]++;
			}else{
				caseCount[this->mGenotypeStr[validCase[_case]][Snp[snp]]]=1;
			}
		}
		for(int _ctrl=0;_ctrl<validCtrl.size();_ctrl++){
			if(ctrlCount.find(this->mGenotypeStr[validCtrl[_ctrl]][Snp[snp]]) != ctrlCount.end()){
				ctrlCount[this->mGenotypeStr[validCtrl[_ctrl]][Snp[snp]]]++;
			}else{
				ctrlCount[this->mGenotypeStr[validCtrl[_ctrl]][Snp[snp]]]=1;
			}
		}

		boost::unordered_map<std::string, int>::iterator iter;
		for(iter=caseCount.begin();iter!=caseCount.end();iter++){
			double rate=(double)iter->second/(double)validCase.size();
			if(0 != rate){
				CaseSum+=rate*log(rate);
			}
//			std::cout<<"case, genotype:"<<iter->first<<",count:"<<iter->second<<",rate="<<rate<<"\n";
		}
		for(iter=ctrlCount.begin();iter!=ctrlCount.end();iter++){
			double rate=(double)iter->second/(double)validCtrl.size();
			if(0 != rate){
				CtrlSum+=rate*log(rate);
			}
//			std::cout<<"ctrl, genotype:"<<iter->first<<",count:"<<iter->second<<",rate="<<rate<<"\n";
		}
	}
	CaseSum*=(-1);
	CtrlSum*=(-1);
}

void GeneInteractionBinary::GenerateGenotypeCombination(std::vector<int> &Snp,std::vector<std::vector<std::string> >& ret){
	std::vector<std::vector<std::string> > idx;
	ret.clear();
	for(int i=0;i<Snp.size();i++){
		std::vector<std::string> _idx;
		boost::unordered_map<std::string, double>::iterator iter;
		for(iter=this->data->vLocusInfo[Snp[i]].BothGenotypeCount.begin();iter!=this->data->vLocusInfo[Snp[i]].BothGenotypeCount.end();iter++){
			_idx.push_back(iter->first);
		}
		idx.push_back(_idx);
	}
	cart_product(ret,idx);
}


} /* namespace SHEsis */
