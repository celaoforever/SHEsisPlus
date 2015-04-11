/*
 * GeneInteractionQTL.cpp
 *
 *  Created on: Jan 10, 2015
 *      Author: ionadmin
 */

#include "GeneInteractionQTL.h"
#include <algorithm>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/stats.hpp>
namespace SHEsis {

typedef boost::accumulators::accumulator_set<double,boost::accumulators::features<boost::accumulators::tag::density> >  acc;
typedef boost::iterator_range<std::vector<std::pair<double,double> >::iterator > histogram_type;

GeneInteractionQTL::GeneInteractionQTL(boost::shared_ptr<SHEsisData> data):GeneInteraction(data)
{}

void GeneInteractionQTL::init(){
	double t=this->getThreshold(this->data->vQuantitativeTrait,100);
//	std::cout<<"Threshold: "<<t<<"\n";
	for(int i=0;i<this->data->getSampleNum();i++){
		if(this->data->vQuantitativeTrait[i]<t)
			this->lowgroup.push_back(i);
		else
			this->highgroup.push_back(i);
	}
//	std::cout<<"sample Num of low group: "<<this->lowgroup.size()<<"\n";
//	std::cout<<"sample Num of high group: "<<this->highgroup.size()<<"\n";
}
GeneInteractionQTL::~GeneInteractionQTL() {
	// TODO Auto-generated destructor stub
}

std::string GeneInteractionQTL::reporthtml(){
	  std::string _res = "\n<h2> Gene Interaction Analysis (QTL): </h2>\n";
	  boost::shared_ptr<SHEsis::CreatHtmlTable> html(new SHEsis::CreatHtmlTable());
	  html->createTable("GeneInteraction_QTL");
	  std::vector<std::string> data;
	  data.push_back("SNP set");
	  data.push_back("Nonmissing");
	  data.push_back("Interaction1");
	  data.push_back("Interaction2");
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
	    data.push_back(convert2string(this->res[i].LowEntropy));
	    data.push_back(convert2string(this->res[i].HighEntropy));
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

std::string GeneInteractionQTL::reporttxt(){
	  std::stringstream _res;
	  _res.precision(3);
	  _res << "\n-------------------------------------------------------\n";
	  _res << "Gene Interaction Analysis (QTL)\n";
	  _res << "-------------------------------------------------------\n";
	  _res << "SNP set\tNonmissing\tdiff\tP value";
	  if (this->adjust) {
		  _res << "\t\tHolm\tSidakSS\tSidakSD\tFDR_BH\tFDR_BY";
	  }
	  _res << "\n";
	  for (int i = 0; i < this->res.size(); i++) {
		  _res << this->res[i].snpset << "\t";
		  _res << this->res[i].nonmissing << "\t";
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


void GeneInteractionQTL::GetNonmissingSamples(std::vector<int>& snp,std::vector<int>& input,std::vector<int>& output ){
	output.clear();
	for(int i=0;i<input.size();i++){
		int sample=input[i];
		bool valid=true;
		for(int i=0;i<snp.size();i++){
			if(this->data->mGenotype[sample][snp[i]][0]==0){
				valid=false;
				break;
			}
		}
		if(valid){
			output.push_back(sample);
		};
	}
}


void GeneInteractionQTL::print(){
//	std::cout<<"snp_set\tnonmissing\tdiff\tpermutation p\tdist p\n";
	for(int i=0;i<this->res.size();i++){
		std::cout<<res[i].snpset<<"\t"<<res[i].diff<</*"\t"<<res[i].p2<<*/"\t"<<res[i].p<<"\n";
	}
}


double GeneInteractionQTL::getThreshold(std::vector<double> vec,int NumOfBins){
	int total=vec.size();
	acc myAccumulator( boost::accumulators::tag::density::num_bins=NumOfBins,
			boost::accumulators::tag::density::cache_size=10);
	for(int j=0;j<vec.size();j++){
		myAccumulator(vec[j]);
	}
	histogram_type hist=boost::accumulators::density(myAccumulator);
	double sum=0;
	for(int i=0;i<hist.size();i++){
		sum+=i*hist[i].second*total;
	}
	double sumB=0;
	double wB=0;
	double wF=0;
	double mB=0;
	double mF=0;
	double max=0;
	double between=0;
	double t1=0;
	double t2=0;
	for(int i=0;i<hist.size();i++){
		wB+=hist[i].second*total;
		if(0 == wB)
			continue;
		wF=total-wB;
		if(0==wF)
			break;
		sumB+=i*hist[i].second*total;
		mB=sumB/wB;
		mF=(sum-sumB)/wF;
		between=wB*wF*(mB-mF)*(mB-mF);
		if(between>=max){
			t1=i;
			if(between>max){
				t2=i;
			}
			max=between;
		}
	}
	return (hist[t1].first+hist[t2].first)/2;

}

void GeneInteractionQTL::CalGeneInteraction(){
	this->init();
	BOOST_ASSERT(this->lowbound<=this->upperbound);
	for(int i=this->lowbound;i<=this->upperbound;i++){
		std::vector<std::vector<int>  > Snp;
		GenerateSNPCombination(i,Snp);
		for(int j=0;j<Snp.size();j++){
			this->res.push_back(this->GetOneSNPCombinationInformationGain(Snp[j]));
		}
	}
//	for(int i=0;i<this->data->getSnpNum();i=i+2){
//		std::vector<int>  Snp;
//		Snp.push_back(i);
//		Snp.push_back(i+1);
////		Snp.push_back(i+2);
//		this->res.push_back(this->GetOneSNPCombinationInformationGain(Snp));
//	}
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

void getStdDevMean(std::vector<double>& v,double & mean,double& stddev){
	mean=0;
	stddev=0;
	for(int i=0;i<v.size();i++){
		mean+=v[i];
	}
	mean/=(double)v.size();
	for(int i=0;i<v.size();i++){
		stddev+=pow(mean-v[i],2.);
	}
	stddev/=(double)(v.size()-1);
	stddev=sqrt(stddev);
}


gxgQTLRes GeneInteractionQTL::GetOneSNPCombinationInformationGain(std::vector<int>& Snp){
	gxgQTLRes ret;
	std::string snpset="";
	for(int i=0;i<Snp.size()-1;i++){
		snpset+=this->data->vLocusName[Snp[i]]+",";
	}
	snpset+=this->data->vLocusName[Snp[Snp.size()-1]];
	ret.snpset=snpset;
	std::cout<<"calculating "<<snpset<<"\n";
	std::vector<int> validSamplesLow;
	std::vector<int> validSamplesHigh;
	this->GetNonmissingSamples(Snp,this->lowgroup,validSamplesLow);
	this->GetNonmissingSamples(Snp,this->highgroup,validSamplesHigh);
	ret.nonmissing=validSamplesLow.size()+validSamplesHigh.size();
	if( validSamplesLow.size()<5 || validSamplesHigh.size()<5){
		std::cout<<"Number of samples is too few, skipping...";
		return ret;
	};
	ret.LowEntropy=this->getInformationInteraction(validSamplesLow,Snp);
	ret.HighEntropy=this->getInformationInteraction(validSamplesHigh,Snp);
	ret.diff=ret.LowEntropy-ret.HighEntropy;
	if(0 == ret.diff){
		ret.p=1;
		return ret;
	}
	if(this->permutation>0){
		std::vector<double> diffs;
		std::vector<int> allSamples=validSamplesLow;
		allSamples.insert(allSamples.end(),validSamplesHigh.begin(),validSamplesHigh.end());
		for(int p=0;p<this->permutation;p++){
			std::random_shuffle(allSamples.begin(),allSamples.end());
			std::vector<int> PermutatedLow(allSamples.begin(),allSamples.begin()+validSamplesLow.size());
			std::vector<int> PermutatedHigh(allSamples.begin()+validSamplesLow.size(),allSamples.end());
			double entropylow=this->getInformationInteraction(PermutatedLow,Snp);
			double entropyhigh=this->getInformationInteraction(PermutatedHigh,Snp);
			diffs.push_back(entropylow-entropyhigh);
		}
//		ret.p2=this->getP(diffs,ret.diff);
		double mean,stddev;
		getStdDevMean(diffs,mean,stddev);
		double t=(ret.diff-mean)/stddev;
		boost::math::normal s(0,1);
		t=t<0?(-1)*t:t;
		try{
			ret.p=2*boost::math::cdf(boost::math::complement(s,t));
		}catch(...){
			ret.p=-999;
		};
	}
	return ret;
}

} /* namespace SHEsis */
