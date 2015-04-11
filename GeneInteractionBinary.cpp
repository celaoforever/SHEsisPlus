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
	  data.push_back("Case Interaction");
	  data.push_back("Control Interaction");
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
	std::cout<<"snp\tctrlII\tcaseII\tdiff\tp permutate\tp dist\tdetail\n";
	for(int i=0;i<this->res.size();i++){
		std::cout<<res[i].snpset<<"\t"<<res[i].ctrlEntropy<<"\t"<<res[i].caseEntropy<<"\t"<<res[i].diff<<"\t"
				/*<<res[i].p2<<"\t"*/<<res[i].p<<"\t";
		std::cout<<"\n";
	}
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


gxgBinaryRes GeneInteractionBinary::GetOneSNPCombinationInformationGain(std::vector<int>& Snp){
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
	ret.caseEntropy =this->getInformationInteraction(validCase,Snp);
	ret.ctrlEntropy=this->getInformationInteraction(validCtrl,Snp);
	ret.diff=ret.caseEntropy-ret.ctrlEntropy;
	ret.nonmissing=validCase.size()+validCtrl.size();
	if(ret.diff==0){
		ret.p=1;
		return ret;
	}
	if(this->permutation>0){
		std::vector<int> allSamples=validCase;
		allSamples.insert(allSamples.end(),validCtrl.begin(),validCtrl.end());
		std::vector<double> diffs;
		for(int p=0;p<this->permutation;p++){
			std::random_shuffle(allSamples.begin(),allSamples.end());
			std::vector<int> PermutatedCase(allSamples.begin(),allSamples.begin()+validCase.size());
			std::vector<int> PermutatedCtrl(allSamples.begin()+validCase.size(),allSamples.end());
			ret.permutatedCaseEntropy.push_back(this->getInformationInteraction(PermutatedCase,Snp));
			ret.permutatedCtrlEntropy.push_back(this->getInformationInteraction(PermutatedCtrl,Snp));
	}
		//stat mean
		for(int i=0;i<this->permutation;i++){
			diffs.push_back(ret.permutatedCaseEntropy[i]-ret.permutatedCtrlEntropy[i]);
			ret.permutatedDiffMean+=(diffs[i]);
		}
		ret.permutatedDiffMean/=(double)this->permutation;
		//stat variance
		for(int i=0;i<this->permutation;i++){
			ret.permutatedDiffVar+=pow((ret.permutatedDiffMean-diffs[i]),2.0);
		}
		ret.permutatedDiffVar/=(double)(this->permutation-1);
//		ret.p2=this->getP(diffs,ret.diff);
		boost::math::normal s(0,1);
		double T=(ret.diff-ret.permutatedDiffMean)/sqrt(ret.permutatedDiffVar);
		T=T<0?(-1)*T:T;
		try{
		ret.p=2*boost::math::cdf(boost::math::complement(s,T));
		}catch(...){
			ret.p=-999;
		}
	}
	return ret;
}

} /* namespace SHEsis */
