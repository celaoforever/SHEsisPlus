/*
 * QTL.cpp
 *
 *  Created on: Oct 3, 2014
 *      Author: ada
 */

#include "QTL.h"
#include "CreatHtmlTable.h"
#include "utility.h"
#include <math.h>
#include <boost/math/distributions/students_t.hpp>
template std::string convert2string<double>(double v);
namespace SHEsis {

QTL::QTL(boost::shared_ptr<SHEsisData> data):data(data),NumOfPermutation(-1) {
	// TODO Auto-generated constructor stub
	this->data->statCount();
}

QTL::~QTL() {
	// TODO Auto-generated destructor stub
}

int QTL::getAlleleCount(int sample,int snp,short allele){
	int count=0;
	for(int p=0;p<this->data->getNumOfChrSet();p++){
		if(this->data->mGenotype[sample][snp][p] == allele)
			count++;
	}
	return count;
}

bool QTL::IsGenotypeMissing(int sample,int snp){
	for(int p=0;p<this->data->getNumOfChrSet();p++){
		if(this->data->mGenotype[sample][snp][p] == GENOTYPE_MISSING)
			return true;
	}
	return false;
}

std::vector<short> QTL::FindAllele(int snp){
	std::vector<short> res;
	boost::unordered_map<short,double>::iterator iter=this->data->vLocusInfo[snp].BothAlleleCount.begin();
	if(this->data->vLocusInfo[snp].BothAlleleCount.size() == 1){
		res.push_back(iter->first);
		return res;
	}
	if(this->data->vLocusInfo[snp].BothAlleleCount.size() == 2){
		double freq=(double)iter->second/(double)this->data->getSampleNum()/(double)this->data->getNumOfChrSet();
		freq<0.5?res.push_back(iter->first):res.push_back((++iter)->first);
		return res;
	};
	for(;iter!=this->data->vLocusInfo[snp].BothAlleleCount.end();++iter){
		res.push_back(iter->first);
	}
	return res;
}
QTLResults QTL::OneLocusQTLAnalysis(int snp,short allele,std::vector<double> qt){
	double g_mean=0,g_var=0;
	double qt_mean=0,qt_var=0;
	double qt_g_covar=0;
	int ValidSampleNum=0;
	for(int iSample=0;iSample<this->data->getSampleNum();iSample++){
		if(0 == qt[iSample]/*this->data->vQuantitativeTrait[iSample]*/)
			continue;
		if(this->IsGenotypeMissing(iSample,snp))
			continue;
		qt_mean+=qt[iSample]/*this->data->vQuantitativeTrait[iSample]*/;
		g_mean+=this->getAlleleCount(iSample,snp,allele);
		ValidSampleNum++;
	}
	qt_mean/=(double)ValidSampleNum;
	g_mean/=(double)ValidSampleNum;

	for(int iSample=0;iSample<this->data->getSampleNum();iSample++){
		if(0 == qt[iSample]/*this->data->vQuantitativeTrait[iSample]*/)
			continue;
		if(this->IsGenotypeMissing(iSample,snp))
			continue;
		qt_var+=(/*this->data->vQuantitativeTrait[iSample]*/qt[iSample]-qt_mean)*
				(/*this->data->vQuantitativeTrait[iSample]*/qt[iSample]-qt_mean);
		double g=this->getAlleleCount(iSample,snp,allele);
		g_var+=(g-g_mean)*(g-g_mean);
		qt_g_covar+=(/*this->data->vQuantitativeTrait[iSample]*/qt[iSample]-qt_mean)*(g-g_mean);
	}

	qt_var/=(double)ValidSampleNum-1;
	g_var/=(double)ValidSampleNum-1;
	qt_g_covar/=(double)ValidSampleNum-1;
	QTLResults res;
	res.beta=qt_g_covar/g_var;
	res.se=sqrt((qt_var/g_var-(qt_g_covar*qt_g_covar)/(g_var*g_var))/((double)ValidSampleNum-2));
	res.T=res.beta/res.se;
	boost::math::students_t dist(ValidSampleNum-2);
	res.p=2*boost::math::cdf(boost::math::complement(dist,fabs(res.T)));
	res.R2=(qt_g_covar*qt_g_covar)/(qt_var*g_var);
	res.allele=allele;
	res.ValidSampleNum=ValidSampleNum;
	return res;
}

double getTheSmallestP(std::vector<QTLResults> res){
	double p=1;
	for(int i=0;i<res.size();i++){
		p=res[i].p<p?res[i].p:p;
	}
	return p;
}


void QTL::QTLPermutation(){
	this->vPermutatedQT=this->data->vQuantitativeTrait;
	for(int i=0;i<this->NumOfPermutation;i++){
	    printf("\rPermutating...%d%%",
	           (int)(100 * (double)i / (double)this->NumOfPermutation));
	    fflush(stdout);
	    std::random_shuffle(this->vPermutatedQT.begin(),
	                        this->vPermutatedQT.end());
	    this->QTLTest(this->vPermutatedQT);
	    double smallestp=getTheSmallestP(this->vResults);
	    this->permutationp.push_back(smallestp);
	}
	printf("\rPermutating...%d%%\n", 100);
	fflush(stdout);

	std::sort(this->permutationp.begin(),this->permutationp.end());
	this->QTLTest();

	for(int i=0;i<this->vResults.size();i++){
		int rank=getRank(this->vResults[i].p,this->permutationp);
		this->vResults[i].permutatedp=(double)rank/(double)this->NumOfPermutation;
	}


}

void QTL::QTLTest(){

	this->vResults.clear();
	for(int iSnp=0;iSnp<this->data->getSnpNum();iSnp++){
		std::vector<short> allelesToTest=this->FindAllele(iSnp);
		QTLResults Res;
		for(int a=0;a<allelesToTest.size();a++){
			QTLResults CurAlleleRes=this->OneLocusQTLAnalysis(iSnp,allelesToTest[a],this->data->vQuantitativeTrait);
			Res=CurAlleleRes<Res?CurAlleleRes:Res;
		};
		this->vResults.push_back(Res);
	}
}

void QTL::QTLTest(std::vector<double> qt){

	this->vResults.clear();
	for(int iSnp=0;iSnp<this->data->getSnpNum();iSnp++){
		std::vector<short> allelesToTest=this->FindAllele(iSnp);
		QTLResults Res;
		for(int a=0;a<allelesToTest.size();a++){
			QTLResults CurAlleleRes=this->OneLocusQTLAnalysis(iSnp,allelesToTest[a],qt);
			Res=CurAlleleRes<Res?CurAlleleRes:Res;
		};
		this->vResults.push_back(Res);
	}
}

std::string QTL::reporthtml() {
  boost::shared_ptr<SHEsis::CreatHtmlTable> html(new SHEsis::CreatHtmlTable());
  html->createTable("QTL_Analysis");
  std::vector<std::string> data;
  data.push_back("SNP");
  data.push_back("Effect allele");
  data.push_back("Nonmissing sample");
  data.push_back("Regression coefficient");
  data.push_back("Standard error");
  data.push_back("Regression r-squared");
  data.push_back("T statistics");
  data.push_back("p value");
  if(this->NumOfPermutation!=-1){
	  data.push_back("permutation p");
  }
  html->addHeadRow(data);
  for (int i = 0; i < this->vResults.size(); i++) {
    data.clear();
    data.push_back(this->data->vLocusName[i]);
    data.push_back(this->data->getallele(this->vResults[i].allele));
    data.push_back(convert2string(this->vResults[i].ValidSampleNum));
    data.push_back(convert2string(this->vResults[i].beta));
    data.push_back(convert2string(this->vResults[i].se));
    data.push_back(convert2string(this->vResults[i].R2));
    data.push_back(convert2string(this->vResults[i].T));
    data.push_back(convert2string(this->vResults[i].p));
    if(this->NumOfPermutation!=-1){
  	  data.push_back(convert2string(this->vResults[i].permutatedp));
    }
    html->addDataRow(data);
  }
  return html->getTable();
}

void QTL::printRes(){
	std::cout<<"\nQTL results:\n";
	for(int i=0;i<this->vResults.size();i++){
		std::cout<<"Locus "<<i<<", nonmissing="<<this->vResults[i].ValidSampleNum<<", allele="<<this->vResults[i].allele
				<<",r2="<<this->vResults[i].R2<<",T="<<this->vResults[i].T<<",beta="<<this->vResults[i].beta<<",se="<<
				this->vResults[i].se<<",p="<<this->vResults[i].p<<",permuation p="<<this->vResults[i].permutatedp<<"\n";
	}
}


} /* namespace SHEsis */
