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
namespace SHEsis {

GeneInteractionQTL::GeneInteractionQTL(boost::shared_ptr<SHEsisData> data):GeneInteraction(data),
		NumBin(2),MinBin(2),MinSamplesPerBin(50){
	this->lr.reset(new linear);

}

GeneInteractionQTL::~GeneInteractionQTL() {
	// TODO Auto-generated destructor stub
}

bool sortSampleByQtl(const qtl2sampleIdx& v1, const qtl2sampleIdx& v2){
	return (v1.qtl<v2.qtl);
}

std::string GeneInteractionQTL::reporthtml(){
	  std::string _res = "\n<h2> Gene Interaction Analysis (QTL): </h2>\n";
	  boost::shared_ptr<SHEsis::CreatHtmlTable> html(new SHEsis::CreatHtmlTable());
	  html->createTable("GeneInteraction_QTL");
	  std::vector<std::string> data;
	  data.push_back("SNP set");
	  data.push_back("Nonmissing");
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


int GeneInteractionQTL::UpdateBins2(std::vector<int>& snp,std::vector<qtl2sampleIdx>& validSamples){
	this->bins.clear();
	std::sort(validSamples.begin(),validSamples.end(),sortSampleByQtl);
	int SamplesPerBin=validSamples.size()/this->NumBin;
	//adjust according to min samples per bin
	int maxbin=validSamples.size()/this->MinSamplesPerBin;
	int curbin=this->NumBin;
	if(curbin>maxbin){
		curbin=maxbin;
		SamplesPerBin=this->MinSamplesPerBin;
	}else{
		int rest=validSamples.size()%curbin;
		if(rest>this->MinSamplesPerBin)
			curbin++;
	}
//	std::cout<<"NumBin="<<curbin<<",samplesPerbin="<<SamplesPerBin<<"\n";
	if(curbin<this->MinBin)
		return 0;

	for(int i=0;i<curbin-1;i++){
		bin b;
		b.meanqtl=0;
		for(int j=0;j<SamplesPerBin;j++){
			int sampleidx=validSamples[i*SamplesPerBin+j].idx;
			double qtl=validSamples[i*SamplesPerBin+j].qtl;
			b.SampleIdx.push_back(sampleidx);
			b.meanqtl+=qtl;
			if(j==0)
				b.start=qtl;
			if(j==SamplesPerBin-1)
				b.end=qtl;

		}
		b.meanqtl/=(double)SamplesPerBin;
		this->bins.push_back(b);
	}

	//deal with the rest
	bin b;
	b.meanqtl=0;
	b.start=validSamples[(curbin-1)*SamplesPerBin].qtl;
	b.end=validSamples[validSamples.size()-1].qtl;
	int c=0;
	for(int j=(curbin-1)*SamplesPerBin;j<validSamples.size();j++){
		b.meanqtl+=validSamples[j].qtl;
		b.SampleIdx.push_back(validSamples[j].idx);
		c++;
	}
	b.meanqtl/=(double)c;
	this->bins.push_back(b);
	return validSamples.size();
}

int GeneInteractionQTL::UpdateBins(std::vector<int>& snp){
	this->bins.clear();
	std::vector<int> validSamples;
	this->GetNonmissingSamples(snp,validSamples);
	//stat mean & stddev
	double mean=0;
	double stddev=0;
	double min=DBL_MAX;
	double max=DBL_MIN;

	for(int i=0;i<validSamples.size();i++)
	{
		mean+=this->data->vQuantitativeTrait[validSamples[i]];
		min=min>this->data->vQuantitativeTrait[validSamples[i]]?this->data->vQuantitativeTrait[validSamples[i]]:min;
		max=max<this->data->vQuantitativeTrait[validSamples[i]]?this->data->vQuantitativeTrait[validSamples[i]]:max;
	}
	mean/=(double)validSamples.size();
	for(int i=0;i<validSamples.size();i++){
		stddev+=pow(mean-this->data->vQuantitativeTrait[validSamples[i]],2.);
	}
	stddev/=(double)(validSamples.size()-1);
	stddev=sqrt(stddev);

	//stat bins
	for(double i=min;i<=max;i+=stddev){
		bin b;
		b.start=i;
		b.end=i+stddev;
		this->statIntervalSampleNum(b,validSamples);
		this->bins.push_back(b);
	}

	//merge bins
	for(std::list<bin>::iterator iter=this->bins.begin();iter!=this->bins.end();iter++){
		while(iter->SampleIdx.size()<this->MinSamplesPerBin){
			if(!this->ForwardMergeBins(iter))
				break;   //reach end of list
		}
	}
	for(std::list<bin>::reverse_iterator iter=this->bins.rbegin();iter!=this->bins.rend();iter++){
		while(iter->SampleIdx.size() < this->MinSamplesPerBin){
			if(!this->BackwardMergeBins(iter))
				break;
		}
	}
	return validSamples.size();
}

void GeneInteractionQTL::GetBins(){
	this->bins.clear();
	//stat mean & stddev
	double mean=0;
	double stddev=0;
	double min=DBL_MAX;
	double max=DBL_MIN;
	for(int i=0;i<this->data->vQuantitativeTrait.size();i++){
		mean+=this->data->vQuantitativeTrait[i];
		min=min>this->data->vQuantitativeTrait[i]?this->data->vQuantitativeTrait[i]:min;
		max=max<this->data->vQuantitativeTrait[i]?this->data->vQuantitativeTrait[i]:max;
	}
	mean/=(double)this->data->vQuantitativeTrait.size();
	for(int i=0;i<this->data->vQuantitativeTrait.size();i++){
		stddev+=pow(mean-this->data->vQuantitativeTrait[i],2.);
	}
	stddev/=(double)(this->data->getSampleNum()-1);
	stddev=sqrt(stddev);

	//stat bins
	for(double i=min;i<=max;i+=stddev){
		bin b;
		b.start=i;
		b.end=i+stddev;
		this->statIntervalSampleNum(b);
		this->bins.push_back(b);
	}

	//merge bins
	for(std::list<bin>::iterator iter=this->bins.begin();iter!=this->bins.end();iter++){
		while(iter->SampleIdx.size()<this->MinSamplesPerBin){
			if(!this->ForwardMergeBins(iter))
				break;   //reach end of list
		}
	}
	for(std::list<bin>::reverse_iterator iter=this->bins.rbegin();iter!=this->bins.rend();iter++){
		while(iter->SampleIdx.size() < this->MinSamplesPerBin){
			if(!this->BackwardMergeBins(iter))
				break;
		}
	}
}

void GeneInteractionQTL::statIntervalSampleNum(bin& b){
	//[start,end)
	for(int i=0;i<this->data->getSampleNum();i++){
		if(this->data->vQuantitativeTrait[i]>=b.start && this->data->vQuantitativeTrait[i]<b.end){
			b.SampleIdx.push_back(i);
		}
	}
}
void GeneInteractionQTL::GetNonmissingSamples2(std::vector<int>& snp,std::vector<qtl2sampleIdx>& output ){
	output.clear();
	for(int sample=0;sample<this->data->getSampleNum();sample++){
		bool valid=true;
		for(int i=0;i<snp.size();i++){
			if(this->data->mGenotype[sample][snp[i]][0]==0){
				valid=false;
				break;
			}
		}
		if(valid){
			qtl2sampleIdx s;
			s.qtl=this->data->vQuantitativeTrait[sample];
			s.idx=sample;
			output.push_back(s);
		};
	}
}
void GeneInteractionQTL::GetNonmissingSamples(std::vector<int>& snp,std::vector<int>& output){
	output.clear();
	for(int sample=0;sample<this->data->getSampleNum();sample++){
		bool valid=true;
		for(int i=0;i<snp.size();i++){
			if(this->data->mGenotype[sample][snp[i]][0]==0){
				valid=false;
				break;
			}
		}
		if(valid)
			output.push_back(sample);
	}
}

void GeneInteractionQTL::statIntervalSampleNum(bin& b,std::vector<int>& samples){
	//[start,end)
	//samples are nonmissing
	for(int i=0;i<samples.size();i++){
		if(this->data->vQuantitativeTrait[samples[i]]>=b.start && this->data->vQuantitativeTrait[samples[i]]<b.end){
			b.SampleIdx.push_back(samples[i]);
		}
	}
}


bool GeneInteractionQTL::ForwardMergeBins(std::list<bin>::iterator& iter){
		//node pointed by iter will be deleted, and merge into iter+1
		if(++iter == this->bins.end())
			return false;
		std::list<bin>::iterator next=iter;
		iter--;
		next->start=iter->start;
		next->SampleIdx.insert(next->SampleIdx.begin(),iter->SampleIdx.begin(),iter->SampleIdx.end());
		this->bins.erase(iter);
		iter=next;
	return true;
}

bool GeneInteractionQTL::BackwardMergeBins(std::list<bin>::reverse_iterator& iter){
	//node pointed by iter will be deleted, and merge into iter-1
	if(++iter == this->bins.rend())
		return false;
	std::list<bin>::reverse_iterator previous=iter;
	iter--;
	previous->end=iter->end;
	previous->SampleIdx.insert(previous->SampleIdx.begin(),iter->SampleIdx.begin(),iter->SampleIdx.end());
	this->bins.erase(--(iter.base()));
	iter=previous;
	return true;
}

void GeneInteractionQTL::print(){
//	std::cout<<"snp_set\tnonmissing\tnum_bin\tpearson correl\tpearson p(permutation)\tpearson p(distribution)\trank correl\trank p(permutation)\trank p(distribution)\trank p (origin)\n";
//	for(int i=0;i<this->res.size();i++){
//		std::cout<<res[i].snpset<<"\t"<<res[i].nonmissing<<"\t"<<res[i].numBin<<"\t"<<res[i].correl<<"\t"<<res[i].p2<<"\t"<<res[i].p<<"\t"<<res[i].rankcorrel<<"\t"<<res[i].rankp2<<"\t"<<res[i].rankp<<"\t"<<res[i].rankporigin<<"\t";
//	       for(int k=0;k<res[i].entropy.size();k++){
//			std::cout<<res[i].entropy[k]<<",";
//		}
//		std::cout<<"\n";
//	}
	std::cout<<"snp_set\tnonmissing\tdiff\tpermutation p\tdist p\n";
	for(int i=0;i<this->res.size();i++){
		std::cout<<res[i].snpset<<"\t"<<res[i].nonmissing<<"\t"<<res[i].diff<<"\t"<<res[i].p2<<"\t"<<res[i].p<<"\n";
	}
}

void GeneInteractionQTL::CalGeneInteraction(){
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

//	for(int i=0;i<this->data->getSnpNum();i=i+2){
//		std::vector<int>  Snp;
//		Snp.push_back(i);
//		Snp.push_back(i+1);
//		this->res.push_back(this->GetOneSNPCombinationInformationGain(Snp));
//	}
}
double GeneInteractionQTL::getIntervalQTL(bin& b){
//	return (b.start+b.end)/2.;
	return b.meanqtl;
}

double normalizeVector(std::vector<double>& v, std::vector<double>& out){
	double mean=0;
	double stddev=0;
	for(int i=0;i<v.size();i++){
		mean+=v[i];
	}
	mean/=(double)v.size();
	for(int i=0;i<v.size();i++){
		stddev+=pow(mean-v[i],2.);
	}
	stddev/=(double)(v.size()-1);
	stddev=sqrt(stddev);
	for(int i=0;i<v.size();i++){
//		double val=(v[i]-mean)/stddev;
		double val=v[i]/2./stddev;
		out.push_back(val);
	}
	return stddev;
}

double getStdDev(std::vector<double>& v,double & mean){
	mean=0;
	double stddev=0;
	for(int i=0;i<v.size();i++){
		mean+=v[i];
	}
	mean/=(double)v.size();
	for(int i=0;i<v.size();i++){
		stddev+=pow(mean-v[i],2.);
	}
	stddev/=(double)(v.size()-1);
	stddev=sqrt(stddev);
	return stddev;
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

void printvector1DD(std::vector<int> snp){
	for(int i=0;i<snp.size();i++){
			std::cout<<snp[i]<<",";
	}
	std::cout<<"\n";
}

void printvector1DD(std::vector<double> snp){
	for(int i=0;i<snp.size();i++){
			std::cout<<snp[i]<<",";
	}
	std::cout<<"\n";
}

void GeneInteractionQTL::randomshuffle(std::vector<qtl2sampleIdx>& samples){
	std::vector<double> qtl;
	for(int i=0;i<samples.size();i++){
		qtl.push_back(samples[i].qtl);
	}
	std::random_shuffle(qtl.begin(),qtl.end());
	for(int i=0;i<samples.size();i++){
		samples[i].qtl=qtl[i];
	}
}

void PairedTtest(std::vector<double>& permutated, std::vector<double>& origin){
	BOOST_ASSERT(permutated.size() == origin.size());
	std::vector<double> t;
	for(int i=0;i<origin.size();i++){
		t.push_back(origin[i]-permutated[i]);
	}
	double mean;
	double stddev=getStdDev(t,mean);
	double T=mean/(stddev/sqrt((double)t.size()));
	int df=t.size()-1;



}

double getEntropyRank(double entropy, std::vector<double> sorted){
	double rankSum=0;
	double count=0;
//	std::cout<<"entropy="<<entropy<<"\n";
	for(int i=0;i<sorted.size();i++){
//		std::cout<<"sorted[i]="<<sorted[i];
		if(entropy==sorted[i]){
			rankSum+=(i+1);
			count++;
//			std::cout<<"==sorted[i]"<<"rankSUm="<<rankSum<<",count="<<count<<"\n";
			continue;
		}
		if(entropy<sorted[i])
			break;
	}
	return (rankSum/count);
}
void GeneInteractionQTL::GetRankCorrelation(std::vector<double> entropy, double& correl,double& p){
	//rank correlation coefficient
	std::vector<double> sortedEntropy=entropy;
	std::vector<double> rank;
	std::sort(sortedEntropy.begin(),sortedEntropy.end());
//	std::cout<<"sorted Entropy:";
//	for(int i=0;i<sortedEntropy.size();i++){
//		std::cout<<sortedEntropy[i]<<",";
//	};
//	std::cout<<"\n";

	for(int i=0;i<entropy.size();i++){
		double rank_=getEntropyRank(entropy[i],sortedEntropy);
		rank.push_back(rank_);
//		std::cout<<"entropy="<<entropy[i]<<",rank="<<rank_<<"\n";
	};
	double d2=0;
	for(int i=0;i<rank.size();i++){
		d2+=(i+1-rank[i])*(i+1-rank[i]);
	}
	double n=(double)rank.size();
	correl=1-6*d2/(n*(n*n-1));
	double t=correl*sqrt((n-2)/(1-correl*correl));
	//student's t distribution with df (n-2)
	try {
		boost::math::students_t dist(n - 2);
		p = 2 * boost::math::cdf(boost::math::complement(dist, fabs(t)));
	}
	catch (...) {
		p = -999;
	}
}
double GeneInteractionQTL::GetCorrelationCoeff(std::vector<double> entropy,std::vector<double> qtl){
	double XVar=0;
	double XMean=0;
	double YVar=0;
	double YMean=0;
	double Covar=0;
	for(int i=0;i<entropy.size();i++){
		XMean+=entropy[i];
	}
	XMean/=(double)entropy.size();
	for(int i=0;i<entropy.size();i++){
			XVar+=pow(entropy[i]-XMean,2);
	}
	XVar=sqrt(XVar);

	for(int i=0;i<qtl.size();i++){
		YMean+=qtl[i];
	}
	YMean/=(double)qtl.size();
	for(int i=0;i<qtl.size();i++){
			YVar+=pow(qtl[i]-YMean,2);
	}
	YVar=sqrt(YVar);
	for(int i=0;i<qtl.size();i++){
		Covar+=(entropy[i]-XMean)*(qtl[i]-YMean);
	}
	double r=Covar/(XVar*YVar);
	return r;

	
}

double GeneInteractionQTL::getP(std::vector<double> permutated,double origin){
	double p=0;
	origin=SHEsisABS(origin);
	for(int i=0;i<permutated.size();i++){
		double abs=SHEsisABS(permutated[i]);
		if(origin<abs){
//			std::cout<<"origin="<<origin<<",permutated="<<abs<<",p++\n";
			p++;
		}
	}
	p=p/(double)permutated.size();

	return p;
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
	std::vector<qtl2sampleIdx> validSamples;
	this->GetNonmissingSamples2(Snp,validSamples);
	ret.nonmissing=this->UpdateBins2(Snp,validSamples);
	if(ret.nonmissing == 0 ){
		std::cout<<"number of bins is less than min bin, skip...\n";
		return ret;
	};
	//calculate information iteraction within each bins
	//std::cout<<"origin\n";
	for(std::list<bin>::iterator iter=this->bins.begin();iter != this->bins.end();iter++){
		ret.entropy.push_back(this->getInformationInteraction(iter->SampleIdx,Snp));
		ret.qtl.push_back(this->getIntervalQTL(*iter));
		//std::cout<<ret.entropy[ret.entropy.size()-1]<<"\t";
	}
	//std::cout<<ret.entropy[0]-ret.entropy[1];
	if(ret.entropy.size()==2){
		double diff=ret.entropy[0]-ret.entropy[1];
		
		if(this->permutation>0){
			//std::cout<<"\npermutating...\n";
			std::vector<double> diffs;
			for(int p=0;p<this->permutation;p++){
				this->randomshuffle(validSamples);
				std::vector<double> entropy;
				std::vector<double> qtl;
				if(this->UpdateBins2(Snp,validSamples)){
					for(std::list<bin>::iterator iter=this->bins.begin();iter != this->bins.end();iter++){
						entropy.push_back(this->getInformationInteraction(iter->SampleIdx,Snp));
						//std::cout<<entropy[entropy.size()-1]<<"\t";
						qtl.push_back(this->getIntervalQTL(*iter));
					}
					//std::cout<<entropy[0]-entropy[1]<<"\n";
				}
				BOOST_ASSERT(entropy.size()==2);
				diffs.push_back(entropy[0]-entropy[1]);
			}
			ret.p2=this->getP(diffs,diff);
			double mean,stddev;
			getStdDevMean(diffs,mean,stddev);
			double t=(diff-mean)/stddev;
			boost::math::normal s(0,1);
			t=t<0?(-1)*t:t;
			try{
			ret.p=2*boost::math::cdf(boost::math::complement(s,t));
			}catch(...){
				ret.p=-999;
			}
			ret.diff=diff;
		}
		return ret;
	}

	//calculate correlation coefficient
	ret.correl=this->GetCorrelationCoeff(ret.entropy,ret.qtl);
	this->GetRankCorrelation(ret.entropy,ret.rankcorrel,ret.rankporigin);
	ret.numBin=this->bins.size();
	if(this->permutation>0){
//		std::cout<<"permutating...\n";
		std::vector<double> correl;
		std::vector<double> rankcorrel;
		for(int p=0;p<this->permutation;p++){
			this->randomshuffle(validSamples);
			std::vector<double> entropy;
			std::vector<double> qtl;
			if(this->UpdateBins2(Snp,validSamples)){
				for(std::list<bin>::iterator iter=this->bins.begin();iter != this->bins.end();iter++){
					entropy.push_back(this->getInformationInteraction(iter->SampleIdx,Snp));
					qtl.push_back(this->getIntervalQTL(*iter));
				}
			}
			double corr=this->GetCorrelationCoeff(entropy,qtl);
			double _corr,_p;
			this->GetRankCorrelation(entropy,_corr,_p);
			rankcorrel.push_back(_corr);
//			std::cout<<ret.snpset<<"\t"<<_corr<<"\n";
			correl.push_back(corr);
		}
		double mean,stddev,rankmean,rankstddev;
		getStdDevMean(correl,mean,stddev);
		getStdDevMean(rankcorrel,rankmean,rankstddev);
		double t=(ret.correl-mean)/stddev;
		double rankt=(ret.rankcorrel-rankmean)/rankstddev;
		boost::math::normal s(0,1);
		t=t<0?(-1)*t:t;
		rankt=rankt<0?(-1)*rankt:rankt;
		ret.p=2*boost::math::cdf(boost::math::complement(s,t));
		ret.p2=this->getP(correl,ret.correl);
		ret.rankp2=this->getP(rankcorrel,ret.rankcorrel);
		ret.rankp=2*boost::math::cdf(boost::math::complement(s,rankt));
//		std::cout<<"mean="<<s.mean()<<",stdev="<<s.standard_deviation()<<",val="<<ret.correl<<",p="<<ret.p<<"\n";
		return ret;
	}
	//linear regression
	this->lr.reset(new linear);
	lr->resetResponse(ret.qtl);
	lr->resetSNP(ret.entropy);
	lr->regress();
	ret.coef=lr->coef(1);
	ret.p=lr->p(1)*2;
	ret.se=lr->se(1);
	return ret;
}




} /* namespace SHEsis */
