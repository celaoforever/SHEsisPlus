/*
 * GeneInteractionBinary.cpp
 *
 *  Created on: Jan 1, 2015
 *      Author: ionadmin
 */

#include "GeneInteractionBinary.h"
#include <algorithm>
#include <boost/math/distributions/normal.hpp>

namespace SHEsis {

GeneInteractionBinary::GeneInteractionBinary(boost::shared_ptr<SHEsisData> data):GeneInteraction(data) {
	// TODO Auto-generated constructor stub
}

GeneInteractionBinary::~GeneInteractionBinary() {
	// TODO Auto-generated destructor stub
}

std::string GeneInteractionBinary::reporthtml(){
	  std::string _res = "\n<h2> Gene Interaction Analysis (Binary): </h2>\n";
	  boost::shared_ptr<SHEsis::CreatHtmlTable> html(new SHEsis::CreatHtmlTable());
	  html->createTable("GeneInteraction_QTL");
	  std::vector<std::string> data;
	  data.push_back("SNP set");
	  data.push_back("Nonmissing");
	  data.push_back("Case Entropy");
	  data.push_back("Control Entropy");
	  data.push_back("diff");
	  data.push_back("p");
	  if (this->adjust) {
	    data.push_back("Holm");
	    data.push_back("SidakSS");
	    data.push_back("SidakSD");
	    data.push_back("FDR_BH");
	    data.push_back("FDR_BY");
	  }
	  html->addHeadRow(data);
	  for (int i = 0; i < this->res.size(); i++) {
	    data.clear();
	    data.push_back(this->res[i].snpset);
	    data.push_back(convert2string(this->res[i].nonmissing));
	    data.push_back(convert2string(this->res[i].caseEntropy));
	    data.push_back(convert2string(this->res[i].ctrlEntropy));
	    data.push_back(convert2string(this->res[i].diff));
	    data.push_back(convert2string(this->res[i].p));
	    if (this->adjust) {
	      data.push_back(convert2string(this->res[i].HolmP));
	      data.push_back(convert2string(this->res[i].SidakSSP));
	      data.push_back(convert2string(this->res[i].SidakSDP));
	      data.push_back(convert2string(this->res[i].BHP));
	      data.push_back(convert2string(this->res[i].BYP));
	    }
	    html->addDataRow(data);
	  }
	  _res += html->getTable();
	  return _res;
}

std::string GeneInteractionBinary::reporttxt(){
	  std::stringstream _res;
	  _res.precision(3);
	  _res << "\n-------------------------------------------------------\n";
	  _res << "Gene Interaction Analysis (Binary)\n";
	  _res << "-------------------------------------------------------\n";
	  _res << "SNP set\tNonmissing\tCase Entropy\tControl Entropy\tDiff\tP value";
	  if (this->adjust) {
		  _res << "\t\tHolm\tSidakSS\tSidakSD\tFDR_BH\tFDR_BY";
	  }
	  _res << "\n";
	  for (int i = 0; i < this->res.size(); i++) {
		  _res << this->res[i].snpset << "\t";
		  _res << this->res[i].nonmissing << "\t";
		  _res << this->res[i].caseEntropy << "\t";
		  _res << this->res[i].ctrlEntropy << "\t";
		  _res << this->res[i].diff << "\t";
		  _res << convert2string(this->res[i].p) << "\t";
	    if (this->adjust) {
	     _res << convert2string(this->res[i].HolmP) << "\t"
	          << convert2string(this->res[i].SidakSSP) << "\t"
	          << convert2string(this->res[i].SidakSDP) << "\t"
	          << convert2string(this->res[i].BHP) << "\t"
	          << convert2string(this->res[i].BYP);
	    }
	    _res << "\n";
	  }
	  _res << "-------------------------------------------------------\n";
	  return _res.str();
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
void GeneInteractionBinary::print(){
	std::cout<<"snp\tctrlII\tcaseII\tdiff\tmean\tvar\tp\tdetail\n";
	for(int i=0;i<this->res.size();i++){
		std::cout<<res[i].snpset<<"\t"<<res[i].ctrlEntropy<<"\t"<<res[i].caseEntropy<<"\t"<<res[i].diff<<"\t"
				<<res[i].permutatedDiffMean<<"\t"<<res[i].permutatedDiffVar<<"\t"<<res[i].p<<"\t";
//		for(int j=0;j<res[i].permutatedCaseEntropy.size();j++){
//				std::cout<<(res[i].permutatedCaseEntropy[j]-res[i].permutatedCtrlEntropy[j])<<",";
//		}
		std::cout<<"\n";
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
			this->res.push_back(this->GetOneSNPCombinationInformationGain2(Snp[j]));
//			this->ThreadPool.enqueue(boost::bind(&GeneInteractionBinary::GetOneSNPCombinationInformationGain,this,Snp[j]));
		}
	}
//		for(int i=0;i<this->data->getSnpNum();i=i+2){
//			std::vector<int>  Snp;
//			Snp.push_back(i);
//			Snp.push_back(i+1);
//			this->res.push_back(this->GetOneSNPCombinationInformationGain2(Snp));
//		}
		 if (this->adjust) {
		    std::vector<MultiComp> originp;
		    std::vector<double> adjusted;
		    for (int i = 0; i < this->res.size(); i++) {
		      MultiComp val;
		      val.p = this->res[i].p;
		      val.idx = i;
		      originp.push_back(val);
		    };
		    std::sort(originp.begin(), originp.end());
		    HolmCorrection(originp, adjusted);
		    for (int i = 0; i < this->res.size(); i++)
		      this->res[i].HolmP = adjusted[i];

		    SidakSDCorrection(originp, adjusted);
		    for (int i = 0; i < this->res.size(); i++)
		      this->res[i].SidakSDP = adjusted[i];

		    SidakSSCorrection(originp, adjusted);
		    for (int i = 0; i < this->res.size(); i++)
		      this->res[i].SidakSSP = adjusted[i];

		    BHCorrection(originp, adjusted);
		    for (int i = 0; i < this->res.size(); i++)
		      this->res[i].BHP = adjusted[i];

		    BYCorrection(originp, adjusted);
		    for (int i = 0; i < this->res.size(); i++)
		      this->res[i].BYP = adjusted[i];
		  }

}


gxgBinaryRes GeneInteractionBinary::GetOneSNPCombinationInformationGain2(std::vector<int>& Snp){
	gxgBinaryRes ret;
	std::string snpset="";
	for(int i=0;i<Snp.size()-1;i++){
		snpset+=this->data->vLocusName[Snp[i]]+",";
	}
	snpset+=this->data->vLocusName[Snp[Snp.size()-1]];
	ret.snpset=snpset;
	std::cout<<"calculating "<<snpset<<"\n";
	std::vector<int> validCase;
	std::vector<int> validCtrl;
	this->getNonmissingSample(Snp,validCase,validCtrl);
//	std::cout<<"nonmissing case:";
//	printvector1D(validCase);
//	std::cout<<"nonmissing ctrl:";
//	printvector1D(validCtrl);
	ret.caseEntropy =this->getInformationInteraction(validCase,Snp);
	ret.ctrlEntropy=this->getInformationInteraction(validCtrl,Snp);
//	std::cout<<"case entropy="<<ret.caseEntropy<<"\n";
//	std::cout<<"ctrl entropy="<<ret.ctrlEntropy<<"\n";
	ret.diff=ret.caseEntropy-ret.ctrlEntropy;
	ret.nonmissing=validCase.size()+validCtrl.size();
	if(ret.diff==0){
		ret.p=1;
		return ret;
	}
//	ret.caseLambda=sqrt(1-exp((-2)*(ret.caseEntropy)));
//	ret.ctrlLambda=sqrt(1-exp((-2)*(ret.ctrlEntropy)));
	if(this->permutation>0){
//		std::cout<<"permutating...\n";
		std::vector<int> allSamples=validCase;
		allSamples.insert(allSamples.end(),validCtrl.begin(),validCtrl.end());
//		std::cout<<"all samples:";
//		printvector1D(allSamples);
		for(int p=0;p<this->permutation;p++){
//			std::cout<<"permutating round "<<p<<"\n";
			std::random_shuffle(allSamples.begin(),allSamples.end());
			std::vector<int> PermutatedCase(allSamples.begin(),allSamples.begin()+validCase.size());
			std::vector<int> PermutatedCtrl(allSamples.begin()+validCase.size(),allSamples.end());
//			std::cout<<"permutated case:";
//			printvector1D(PermutatedCase);
//			std::cout<<"permutated ctrl:";
//			printvector1D(PermutatedCtrl);
			ret.permutatedCaseEntropy.push_back(this->getInformationInteraction(PermutatedCase,Snp));
			ret.permutatedCtrlEntropy.push_back(this->getInformationInteraction(PermutatedCtrl,Snp));
//			std::cout<<"permutated case entropy:"<<ret.permutatedCaseEntropy[ret.permutatedCaseEntropy.size()-1]<<"\n";
//			std::cout<<"permutated ctrl entropy:"<<ret.permutatedCtrlEntropy[ret.permutatedCtrlEntropy.size()-1]<<"\n";
		}
		//stat mean
		for(int i=0;i<this->permutation;i++){
			ret.permutatedDiffMean+=(ret.permutatedCaseEntropy[i]-ret.permutatedCtrlEntropy[i]);
		}
		ret.permutatedDiffMean/=(double)this->permutation;
		//stat variance
		for(int i=0;i<this->permutation;i++){
			ret.permutatedDiffVar+=pow((ret.permutatedDiffMean-(ret.permutatedCaseEntropy[i]-ret.permutatedCtrlEntropy[i])),2.0);
		}
		ret.permutatedDiffVar/=(double)(this->permutation-1);
		boost::math::normal s(0,1);
		double T=(ret.diff-ret.permutatedDiffMean)/sqrt(ret.permutatedDiffVar);
		T=T<0?(-1)*T:T;
		try{
		ret.p=boost::math::cdf(boost::math::complement(s,T));
		}catch(...){
			ret.p=-999;
		}
//		std::cout<<"mean="<<s.mean()<<",stdev="<<s.standard_deviation()<<",abs_val="<<_abs<<",p="<<ret.p<<"\n";
	}
	return ret;
}

//not used
gxgBinaryRes GeneInteractionBinary::GetOneSNPCombinationInformationGain(std::vector<int>& Snp){
	gxgBinaryRes ret;
	std::string snpset="";
	for(int i=0;i<Snp.size()-1;i++){
		snpset+=this->data->vLocusName[Snp[i]]+",";
	}
	snpset+=this->data->vLocusName[Snp[Snp.size()-1]];
	ret.snpset=snpset;
	std::cout<<"calculating "<<snpset<<"\n";
	double caseEntropy=0,controlEntropy=0;
	std::vector<std::vector<std::string> > cp;
	this->GenerateGenotypeCombination(Snp,cp);
	this->GetInformationGain(Snp,cp,ret.caseEntropy,ret.ctrlEntropy);
	ret.diff=ret.caseEntropy-ret.ctrlEntropy;
//	ret.caseLambda=sqrt(1-exp((-2)*(ret.caseEntropy)));
//	ret.ctrlLambda=sqrt(1-exp((-2)*(ret.ctrlEntropy)));
	return ret;
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
	std::cout<<"case Entropy sum="<<caseSum<<"\n";
	std::cout<<"ctrl Entropy sum="<<ctrlSum<<"\n";
	//calculate case information gain
	double caseMutual=0;
	double ctrlMutual=0;
	for(int cpIdx=0;cpIdx<cp.size();cpIdx++){
		int count=0;
		for(int sample=0;sample<validCase.size();sample++){
			bool equal=true;
			for(int i=0;i<cp[cpIdx].size();i++){
				if(!this->genotypeEqual(cp[cpIdx][i],validCase[sample],Snp[i])){
					equal = false;
					break;
				}
			}
			if(equal)
				count++;
		}
		std::cout<<"genotype combination for case:\n";
		for(int i=0;i<cp[cpIdx].size();i++){
			std::cout<<cp[cpIdx][i]<<",";
		}

		double rate=(double)count/(double)validCase.size();
		std::cout<<"\nrate="<<rate<<"\n";
		if(rate!=0)
			caseMutual+=rate*log(rate);
	}

	//calculate case information gain
	for(int cpIdx=0;cpIdx<cp.size();cpIdx++){
		int count=0;
		for(int sample=0;sample<validCtrl.size();sample++){
			bool equal=true;
			for(int i=0;i<cp[cpIdx].size();i++){
				if(! this->genotypeEqual(cp[cpIdx][i],validCtrl[sample],Snp[i])){
					equal = false;
					break;
				}
			}
			if(equal)
				count++;
		}
		std::cout<<"genotype combination for control:\n";
		for(int i=0;i<cp[cpIdx].size();i++){
			std::cout<<cp[cpIdx][i]<<",";
		}

		double rate=(double)count/(double)validCtrl.size();
		std::cout<<"\nrate="<<rate<<"\n";
		if(rate!=0)
			ctrlMutual+=rate*log(rate);
	}
	caseMutual*=(-1);
	ctrlMutual*=(-1);
	std::cout<<"case Mutual entropy="<<caseMutual<<"\n";
	std::cout<<"ctrl Mutual entropy="<<ctrlMutual<<"\n";
	caseGain=caseSum-caseMutual;
	ctrlGain=ctrlSum-ctrlMutual;
}




//not used
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
			std::cout<<"case, genotype:"<<iter->first<<",count:"<<iter->second<<",rate="<<rate<<"\n";
		}
		for(iter=ctrlCount.begin();iter!=ctrlCount.end();iter++){
			double rate=(double)iter->second/(double)validCtrl.size();
			if(0 != rate){
				CtrlSum+=rate*log(rate);
			}
			std::cout<<"ctrl, genotype:"<<iter->first<<",count:"<<iter->second<<",rate="<<rate<<"\n";
		}
	}
	CaseSum*=(-1);
	CtrlSum*=(-1);
}


} /* namespace SHEsis */
