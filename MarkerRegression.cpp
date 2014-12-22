/*
 * MarkerRegression.cpp
 *
 *  Created on: Dec 19, 2014
 *      Author: ionadmin
 */

#include "MarkerRegression.h"
#include "utility.h"
namespace SHEsis {
MarkerRegression::MarkerRegression(boost::shared_ptr<SHEsisData> d):data(d),model(ADDICTIVE),adjust(false),permutation(-1)
{
	// TODO Auto-generated constructor stub


}

MarkerRegression::~MarkerRegression() {
	// TODO Auto-generated destructor stub
}
double MarkerRegression::codeAllele(int sample, int snp, short allele){
	int count=this->getAlleleCount(sample,snp,allele);
	double code;
	switch(this->model){
	case ADDICTIVE:
		code=count;
		break;
	case DOMINANT:
		code=count>0?1:0;
		break;
	case RECESSIVE:
		code=count==this->data->getNumOfChrSet()?1:0;
		break;
	default:
		BOOST_ASSERT(1==0);
		break;
	}
	return code;
}

double getTheSmallestP(std::vector<RegressionRes> res) {
  double p = 1;
  for (int i = 0; i < res.size(); i++) {
	  if(res[i].p == -999)
		  continue;
    p = res[i].p < p ? res[i].p : p;
  }
  return p;
}
void MarkerRegression::RegressPermutation(std::vector<double>& v){
	std::vector<double> permutatedp;
	  for (int i = 0; i < this->permutation; i++) {
	    printf("\rPermutating...%d%%",
	           (int)(100 * (double)i / (double)this->permutation));
	    fflush(stdout);
	    std::random_shuffle(v.begin(), v.end());
	    this->lr->setReponse(v);
		std::vector<RegressionRes> PermutateRes;
		for(int i=0;i<this->data->getSnpNum();i++){
			 std::vector<short> allelesToTest = this->FindAllele(i);
			 RegressionRes Res;
			 if(allelesToTest.size()==0){
				 Res.p=-999;
				 PermutateRes.push_back(Res);
				 continue;
			 }
		      for (int a = 0; a < allelesToTest.size(); a++) {
		    	  RegressionRes CurAlleleRes = this->OneLocusRegression(i, allelesToTest[a]);
		        Res = CurAlleleRes < Res ? CurAlleleRes : Res;
		      };
		      PermutateRes.push_back(Res);
		}

	    double smallestp = getTheSmallestP(PermutateRes);
	    permutatedp.push_back(smallestp);
	  }
	  printf("\rPermutating...%d%%\n", 100);
	  fflush(stdout);

	  std::sort(permutatedp.begin(), permutatedp.end());
	  for (int i = 0; i < this->vResults.size(); i++) {
	    int rank = getRank(this->vResults[i].p, permutatedp);
	    this->vResults[i].permutationP =
	        this->vResults[i].p > 0 ? (double)rank / (double)this->permutation
	                                : -999;
	  }

}

void MarkerRegression::regressAll(){
	//reset regression handler
	std::vector<double> resp;
	if(this->data->vLabel.size()==0 && this->data->vQuantitativeTrait.size()>0){
		this->lr.reset(new linear);
		resp=this->data->vQuantitativeTrait;
//		this->lr->setReponse(this->data->vQuantitativeTrait);
		this->data->statCount();
	}else if (this->data->vLabel.size()>0 && this->data->vQuantitativeTrait.size()==0){
		this->lr.reset(new logistic);
		this->data->statCount(this->data->vLabel);
//		std::vector<double> binary;
		for(int i=0;i<this->data->vLabel.size();i++){
			if(this->data->vLabel[i] == CASE)
				resp.push_back(2);
			else
				resp.push_back(1);
		}
//		this->lr->setReponse(binary);
	}else
		BOOST_ASSERT(1==0);
	this->lr->setReponse(resp);
	//set covar
	if(this->data->covar.size()!=0)
		lr->setCovar(this->data->covar);

	this->vResults.clear();
	for(int i=0;i<this->data->getSnpNum();i++){
		 std::vector<short> allelesToTest = this->FindAllele(i);
		 RegressionRes Res;
		 if(allelesToTest.size()==0){
			 Res.p=-999;
			 this->vResults.push_back(Res);
			 continue;
		 }
	      for (int a = 0; a < allelesToTest.size(); a++) {
	    	  RegressionRes CurAlleleRes = this->OneLocusRegression(i, allelesToTest[a]);
	        Res = CurAlleleRes < Res ? CurAlleleRes : Res;
	      };
	      this->vResults.push_back(Res);
	}
	if(this->permutation!=-1){
		this->RegressPermutation(resp);
	}

	 if (this->adjust) {
	    std::vector<MultiComp> originp;
	    std::vector<double> adjusted;
	    for (int i = 0; i < this->vResults.size(); i++) {
	      MultiComp val;
	      val.p = this->vResults[i].p;
	      val.idx = i;
	      originp.push_back(val);
	    };
	    std::sort(originp.begin(), originp.end());
	    HolmCorrection(originp, adjusted);
	    for (int i = 0; i < this->vResults.size(); i++)
	      this->vResults[i].HolmP = adjusted[i];

	    SidakSDCorrection(originp, adjusted);
	    for (int i = 0; i < this->vResults.size(); i++)
	      this->vResults[i].SidakSDP = adjusted[i];

	    SidakSSCorrection(originp, adjusted);
	    for (int i = 0; i < this->vResults.size(); i++)
	      this->vResults[i].SidakSSP = adjusted[i];

	    BHCorrection(originp, adjusted);
	    for (int i = 0; i < this->vResults.size(); i++)
	      this->vResults[i].BHP = adjusted[i];

	    BYCorrection(originp, adjusted);
	    for (int i = 0; i < this->vResults.size(); i++)
	      this->vResults[i].BYP = adjusted[i];
	  }
}

RegressionRes MarkerRegression::OneLocusRegression(int snp,short allele){
	double ave=0;
	int nonmissing=0;
	//replace missing GENOTYPE with average;
	std::vector<double> SNP(this->data->getSampleNum());
	for(int sample=0;sample<this->data->getSampleNum();sample++){
		if(!this->IsGenotypeMissing(sample,snp)){
			double code=this->codeAllele(sample,snp,allele);
			SNP[sample]=code;
			ave+=code;
			nonmissing++;
		}else{
			SNP[sample]=-1;
		}
	}
	ave/=(double)nonmissing;
	for(int i=0;i<SNP.size();i++){
			if(SNP[i] == -1)
				SNP[i]=ave;
	}
	//set snp
	lr->resetSNP(SNP);
	lr->regress();
	RegressionRes r;
	r.allele=allele;
	r.coef=lr->coef(1);
	r.p=lr->p(1);
	r.se=lr->se(1);
	r.nonmissing=nonmissing;
	return r;
}

int MarkerRegression::getAlleleCount(int sample, int snp, short allele) {
  int count = 0;
  for (int p = 0; p < this->data->getNumOfChrSet(); p++) {
    if (this->data->mGenotype[sample][snp][p] == allele) count++;
  }
  return count;
}

bool MarkerRegression::IsGenotypeMissing(int sample, int snp) {
  for (int p = 0; p < this->data->getNumOfChrSet(); p++) {
    if (this->data->mGenotype[sample][snp][p] == GENOTYPE_MISSING) return true;
  }
  return false;
}

std::vector<short> MarkerRegression::FindAllele(int snp) {
  std::vector<short> res;
  boost::unordered_map<short, double>::iterator iter =
      this->data->vLocusInfo[snp].BothAlleleCount.begin();
  if (this->data->vLocusInfo[snp].BothAlleleCount.size() == 1) {
    //		res.push_back(iter->first);
    return res;
  }
  if (this->data->vLocusInfo[snp].BothAlleleCount.size() == 2) {
    double freq = (double)iter->second / (double)this->data->getSampleNum() /
                  (double)this->data->getNumOfChrSet();
    freq < 0.5 ? res.push_back(iter->first) : res.push_back((++iter)->first);
    return res;
  };
  for (; iter != this->data->vLocusInfo[snp].BothAlleleCount.end(); ++iter) {
    res.push_back(iter->first);
  }
  return res;
}

std::string MarkerRegression::reporthtml(){
	  std::string res = "\n<h2> Single Locus Association Test (Regression): </h2>\n";
	  boost::shared_ptr<SHEsis::CreatHtmlTable> html(new SHEsis::CreatHtmlTable());
	  html->createTable("Marker_Regression");
	  std::vector<std::string> data;
	  data.push_back("SNP");
	  data.push_back("Effect allele");
	  data.push_back("Nonmissing");
	  data.push_back("coeff");
	  data.push_back("SE");
	  data.push_back("p");
	  if(this->permutation!=-1){
		  data.push_back("Permutation P");
	  }
	  if (this->adjust) {
	    data.push_back("Holm");
	    data.push_back("SidakSS");
	    data.push_back("SidakSD");
	    data.push_back("FDR_BH");
	    data.push_back("FDR_BY");
	  }
	  html->addHeadRow(data);
	  for (int i = 0; i < this->vResults.size(); i++) {
	    data.clear();
	    data.push_back(this->data->vLocusName[i]);
	    data.push_back(this->data->getallele(this->vResults[i].allele));
	    data.push_back(convert2string(this->vResults[i].nonmissing));
	    data.push_back(convert2string(this->vResults[i].coef));
	    data.push_back(convert2string(this->vResults[i].se));
	    data.push_back(convert2string(this->vResults[i].p));
		if(this->permutation!=-1){
			  data.push_back(convert2string(this->vResults[i].permutationP));
		}
	    if (this->adjust) {
	      data.push_back(convert2string(this->vResults[i].HolmP));
	      data.push_back(convert2string(this->vResults[i].SidakSSP));
	      data.push_back(convert2string(this->vResults[i].SidakSDP));
	      data.push_back(convert2string(this->vResults[i].BHP));
	      data.push_back(convert2string(this->vResults[i].BYP));
	    }
	    html->addDataRow(data);
	  }
	  res += html->getTable();
	  return res;
}

std::string MarkerRegression::reporttxt() {
  std::stringstream res;
  res.precision(3);
  res << "\n-------------------------------------------------------\n";
  res << "Single Locus Association Test (Regression)\n";
  res << "-------------------------------------------------------\n";
  res << "SNP\tEffect allele\tNonmissing\tRegression coefficient\t\tStandard "
         "error\tP value";
  if (this->adjust) {
    res << "\t\tHolm\tSidakSS\tSidakSD\tFDR_BH\tFDR_BY";
  }
  res << "\n";
  for (int i = 0; i < this->vResults.size(); i++) {
    res << this->data->vLocusName[i] << "\t";
    res << this->data->getallele(this->vResults[i].allele) << "\t\t";
    res << this->vResults[i].nonmissing << "\t\t";
    res << this->vResults[i].coef << "\t\t\t";
    res << this->vResults[i].se << "\t\t";
    res << convert2string(this->vResults[i].p) << "\t\t";
    if (this->adjust) {
      res << convert2string(this->vResults[i].HolmP) << "\t"
          << convert2string(this->vResults[i].SidakSSP) << "\t"
          << convert2string(this->vResults[i].SidakSDP) << "\t"
          << convert2string(this->vResults[i].BHP) << "\t"
          << convert2string(this->vResults[i].BYP);
    }
    res << "\n";
  }
  res << "-------------------------------------------------------\n";
  return res.str();
}


}
