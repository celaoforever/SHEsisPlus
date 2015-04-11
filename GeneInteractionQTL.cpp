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

GeneInteractionQTL::GeneInteractionQTL(boost::shared_ptr<SHEsisData> data):GeneInteraction(data)
{

}
void GeneInteractionQTL::init(){
//	std::vector<double> sortedQTL=this->data->vQuantitativeTrait;
//	std::sort(sortedQTL.begin(),sortedQTL.end());
//	double medium=0;
//	if(0==sortedQTL.size()%2){
//		 int idx1=sortedQTL.size()/2-1;
//		 int idx2=sortedQTL.size()/2;
//		 medium=(sortedQTL[idx1]+sortedQTL[idx2])/2;
//	}else{
//		int idx=(sortedQTL.size()-1)/2;
//		medium=sortedQTL[idx];
//	}
//	std::vector<int> equalMedium;
//	for(int i=0;i<this->data->getSampleNum();i++){
//		if(this->data->vQuantitativeTrait[i]<medium){
//			this->lowgroup.push_back(i);
//		}else if(this->data->vQuantitativeTrait[i]>medium){
//			this->highgroup.push_back(i);
//		}else{
//			equalMedium.push_back(i);
//		}
//	}
//	if(this->lowgroup.size()<this->highgroup.size()){
//		this->lowgroup.insert(this->lowgroup.end(),equalMedium.begin(),equalMedium.end());
//	}else{
//		this->highgroup.insert(this->highgroup.end(),equalMedium.begin(),equalMedium.end());
//	}
	double t=this->getThreshold(this->data->vQuantitativeTrait,100);
	std::cout<<"Threshold: "<<t<<"\n";
	for(int i=0;i<this->data->getSampleNum();i++){
		if(this->data->vQuantitativeTrait[i]<t)
			this->lowgroup.push_back(i);
		else
			this->highgroup.push_back(i);
	}
	std::cout<<"sample Num of low group: "<<this->lowgroup.size()<<"\n";
	std::cout<<"sample Num of high group: "<<this->highgroup.size()<<"\n";

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
		std::cout<<res[i].snpset<<"\t"<<res[i].diff<<"\t"<<res[i].p2<<"\t"<<res[i].p<<"\n";
	}
}

typedef boost::accumulators::accumulator_set<double,boost::accumulators::features<boost::accumulators::tag::density> >  acc;
typedef boost::iterator_range<std::vector<std::pair<double,double> >::iterator > histogram_type;
void printvector1D(std::vector<double> snp){
	for(int i=0;i<snp.size();i++){
			std::cout<<snp[i]<<",";
	}
	std::cout<<"\n";
}

double GeneInteractionQTL::getThreshold(std::vector<double> vec,int NumOfBins){
//	printvector1D(this->data->vQuantitativeTrait);
	int total=vec.size();
	acc myAccumulator( boost::accumulators::tag::density::num_bins=NumOfBins,
			boost::accumulators::tag::density::cache_size=10);
	for(int j=0;j<vec.size();j++){
		myAccumulator(vec[j]);
	}
	histogram_type hist=boost::accumulators::density(myAccumulator);
//	for(int i=0;i<hist.size();i++){
//		std::cout<<"Bin lower bound: "<<hist[i].first<<",value: "<<hist[i].second<<"\n";
//	}
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

bool sortSampleByQtl(const qtl2sampleIdx& v1, const qtl2sampleIdx& v2){
return (v1.qtl<v2.qtl);
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
// std::cout<<"NumBin="<<curbin<<",samplesPerbin="<<SamplesPerBin<<"\n";
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

double GeneInteractionQTL::getIntervalQTL(bin& b){
// return (b.start+b.end)/2.;
return b.meanqtl;
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

double getEntropyRank(double entropy, std::vector<double> sorted){
double rankSum=0;
double count=0;
// std::cout<<"entropy="<<entropy<<"\n";
for(int i=0;i<sorted.size();i++){
// std::cout<<"sorted[i]="<<sorted[i];
if(entropy==sorted[i]){
rankSum+=(i+1);
count++;
// std::cout<<"==sorted[i]"<<"rankSUm="<<rankSum<<",count="<<count<<"\n";
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
// std::cout<<"sorted Entropy:";
// for(int i=0;i<sortedEntropy.size();i++){
// std::cout<<sortedEntropy[i]<<",";
// };
// std::cout<<"\n";
for(int i=0;i<entropy.size();i++){
double rank_=getEntropyRank(entropy[i],sortedEntropy);
rank.push_back(rank_);
// std::cout<<"entropy="<<entropy[i]<<",rank="<<rank_<<"\n";
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

double getStdDev(std::vector<double>& v){
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
return stddev;
}

double getBinStdDev(std::list<bin> &_bins){
	double mean=0;
	double stddev=0;
	for(std::list<bin>::iterator iter=_bins.begin();iter!=_bins.end();iter++){
	mean+=iter->entropy;
	}
	mean/=(double)_bins.size();
	for(std::list<bin>::iterator iter=_bins.begin();iter!=_bins.end();iter++){
	stddev+=pow(mean-iter->entropy,2.);
	}
	stddev/=(double)(_bins.size()-1);
	stddev=sqrt(stddev);
	return stddev;
}

bool GeneInteractionQTL::ForwardMergeBins(std::list<bin>::iterator& iter,std::vector<int>& Snps){
//node pointed by iter will be deleted, and merge into iter+1
if(++iter == this->bins.end())
return false;
std::list<bin>::iterator next=iter;
iter--;
next->start=iter->start;
next->SampleIdx.insert(next->SampleIdx.begin(),iter->SampleIdx.begin(),iter->SampleIdx.end());
this->bins.erase(iter);
iter=next;
iter->entropy=this->getInformationInteraction(iter->SampleIdx,Snps);
return true;
}

void printlist(std::list<bin> _bins){
	for(std::list<bin>::iterator iter=_bins.begin();iter!=_bins.end();iter++){
		std::cout<<"sample num:"<<iter->SampleIdx.size()<<",entropy:"<<iter->entropy<<"\n";
	}
	std::cout<<"\n";
}
std::list<bin>::iterator GeneInteractionQTL::findMin(std::list<bin>& _bins,double& mindiff){
	std::list<bin>::iterator ret=_bins.end();
	mindiff=1;
	for(std::list<bin>::iterator iter=_bins.begin();iter!=_bins.end();iter++){
		if(++iter==_bins.end())
			break;
		double nextEn=(iter)->entropy;
		double thisEn=(--iter)->entropy;
		double curdiff=SHEsisABS(nextEn-thisEn);
		if(curdiff<mindiff){
			ret=iter;
			mindiff=curdiff;
		}
	}
	return ret;
}

void GeneInteractionQTL::MergeBins(std::vector<int>& Snp,std::list<bin>& _bins,int count,std::vector<double>& entropy){

	if(count==0)
		return;
	entropy.clear();
	std::vector<double> qtl;
	for(std::list<bin>::iterator iter=_bins.begin();iter != _bins.end();iter++){
		iter->entropy=this->getInformationInteraction(iter->SampleIdx,Snp);
		entropy.push_back(iter->entropy);
		qtl.push_back(iter->meanqtl);
	}
	this->lr.reset(new linear);
	lr->resetSNP(entropy);
	lr->resetResponse(qtl);
	lr->regress();
	double coef=lr->coef(1);
	double p=lr->p(1)*2;
	double se=lr->se(1);
	std::cout<<"coef="<<coef<<",se="<<se<<",p="<<p<<"\n";
	double diff;
	double threshold=.5*getStdDev(entropy);
	std::cout<<"origin:\n";
	printlist(_bins);
	if(p>0.5)
		return;
	std::cout<<"merging start:\n";
	for(std::list<bin>::iterator iter=this->findMin(_bins,diff);iter!=_bins.end()&&diff<threshold;iter=this->findMin(_bins,diff))
	{
		if(!this->ForwardMergeBins(iter,Snp))
			break; //reach end of list
//		threshold=.5*getBinStdDev(_bins);
		printlist(_bins);
	}
	entropy.clear();
	for(std::list<bin>::iterator iter=_bins.begin();iter != _bins.end();iter++){
		iter->entropy=this->getInformationInteraction(iter->SampleIdx,Snp);
		entropy.push_back(iter->entropy);
	}
	std::cout<<"merging end:\n";
	return;

	std::cout<<"before merge: ";
	for(std::list<bin>::iterator iter=_bins.begin();iter!=_bins.end();iter++){
		std::cout<<iter->entropy<<",";
	}
	std::cout<<"\n";
	double stddev=getStdDev(entropy);
	int numStd=1;
	count=0;
	std::cout<<"merge start\n";
	for(std::list<bin>::iterator iter=_bins.begin();iter!=_bins.end();iter++){
		double thisEntropy=(iter)->entropy;
		if(++iter==_bins.end())
			break;
		double nextEntropy=(iter)->entropy;
		iter--;
		printlist(_bins);
		std::cout<<"\n";
//		std::cout<<"this:"<<thisEntropy<<",next:"<<nextEntropy<<",threshold:"<<numStd*stddev<<"\n";
		while((SHEsisABS((thisEntropy-nextEntropy)))<numStd*stddev){
			std::cout<<"this:"<<thisEntropy<<",next:"<<nextEntropy<<"diff:"<<(SHEsisABS((thisEntropy-nextEntropy)))<<",threshold:"<<numStd*stddev<<"\n";
			if(!this->ForwardMergeBins(iter,Snp))
				break; //reach end of list
			count++;
			printlist(_bins);
			std::cout<<"\n";
			thisEntropy=(iter)->entropy;
			if(++iter==_bins.end()){
				iter--;
				break;
			}
			nextEntropy=(iter)->entropy;
			iter--;


		}
	}
	std::cout<<"merge end\n";
	this->MergeBins(Snp,_bins,count,entropy);




}
gxgQTLRes GeneInteractionQTL::GetOneSNPCombinationInformationGain2(std::vector<int>& Snp){
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
//for(std::list<bin>::iterator iter=this->bins.begin();iter != this->bins.end();iter++){
//ret.entropy.push_back(this->getInformationInteraction(iter->SampleIdx,Snp));
//ret.qtl.push_back(this->getIntervalQTL(*iter));
////std::cout<<ret.entropy[ret.entropy.size()-1]<<"\t";
//}
this->MergeBins(Snp,this->bins,1,ret.entropy);
//std::cout<<ret.entropy[0]-ret.entropy[1];
//if(ret.entropy.size()==2){
//double diff=ret.entropy[0]-ret.entropy[1];
//if(this->permutation>0){
////std::cout<<"\npermutating...\n";
//std::vector<double> diffs;
//for(int p=0;p<this->permutation;p++){
//this->randomshuffle(validSamples);
//std::vector<double> entropy;
//std::vector<double> qtl;
//if(this->UpdateBins2(Snp,validSamples)){
//for(std::list<bin>::iterator iter=this->bins.begin();iter != this->bins.end();iter++){
//entropy.push_back(this->getInformationInteraction(iter->SampleIdx,Snp));
////std::cout<<entropy[entropy.size()-1]<<"\t";
//qtl.push_back(this->getIntervalQTL(*iter));
//}
////std::cout<<entropy[0]-entropy[1]<<"\n";
//}
//BOOST_ASSERT(entropy.size()==2);
//diffs.push_back(entropy[0]-entropy[1]);
//}
//ret.p2=this->getP(diffs,diff);
//double mean,stddev;
//getStdDevMean(diffs,mean,stddev);
//double t=(diff-mean)/stddev;
//boost::math::normal s(0,1);
//t=t<0?(-1)*t:t;
//try{
//ret.p=2*boost::math::cdf(boost::math::complement(s,t));
//}catch(...){
//ret.p=-999;
//}
//ret.diff=diff;
//}
//return ret;
//}
//calculate correlation coefficient
this->GetRankCorrelation(ret.entropy,ret.rankcorrel,ret.rankporigin);
std::cout<<"rank correl:"<<ret.rankcorrel<<",rank p:"<<ret.rankporigin<<"\n";
ret.numBin=this->bins.size();
if(this->permutation>0){
 std::cout<<"permutating...\n";
std::vector<double> correl;
std::vector<double> rankcorrel;
for(int p=0;p<this->permutation;p++){
std::cout<<"round "<<p<<"\n";
this->randomshuffle(validSamples);
std::vector<double> entropy;
std::vector<double> qtl;
if(this->UpdateBins2(Snp,validSamples)){
//for(std::list<bin>::iterator iter=this->bins.begin();iter != this->bins.end();iter++){
//entropy.push_back(this->getInformationInteraction(iter->SampleIdx,Snp));
//qtl.push_back(this->getIntervalQTL(*iter));
//}
	this->MergeBins(Snp,this->bins,1,entropy);
}
double _corr,_p;
this->GetRankCorrelation(entropy,_corr,_p);
std::cout<<"rank correl:"<<_corr<<",rank p:"<<_p<<"\n";
rankcorrel.push_back(_corr);
// std::cout<<ret.snpset<<"\t"<<_corr<<"\n";
}
double rankmean,rankstddev;
getStdDevMean(rankcorrel,rankmean,rankstddev);
double rankt=(ret.rankcorrel-rankmean)/rankstddev;
boost::math::normal s(0,1);
rankt=rankt<0?(-1)*rankt:rankt;
ret.rankp2=this->getP(rankcorrel,ret.rankcorrel);
try{
ret.rankp=2*boost::math::cdf(boost::math::complement(s,rankt));
}catch(...){
	ret.rankp=-999;
}
std::cout<<"rank t="<<rankt<<"\n";
std::cout<<"permutate_p="<<ret.rankp2<<",dist_p="<<ret.rankp<<"\n";
return ret;
}
return ret;
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
//	std::cout<<"origin:"<<ret.diff<<"\n";
	if(0 == ret.diff){
		ret.p=1;
		ret.p2=1;
		return ret;
	}
	if(this->permutation>0){
//		std::cout<<"permutating...\n";
		std::vector<double> diffs;
		std::vector<int> allSamples=validSamplesLow;
		allSamples.insert(allSamples.end(),validSamplesHigh.begin(),validSamplesHigh.end());
		for(int p=0;p<this->permutation;p++){
			std::random_shuffle(allSamples.begin(),allSamples.end());
			std::vector<int> PermutatedLow(allSamples.begin(),allSamples.begin()+validSamplesLow.size());
			std::vector<int> PermutatedHigh(allSamples.begin()+validSamplesLow.size(),allSamples.end());
			double entropylow=this->getInformationInteraction(PermutatedLow,Snp);
			double entropyhigh=this->getInformationInteraction(PermutatedHigh,Snp);
//			std::cout<<(entropylow-entropyhigh)<<",";
			diffs.push_back(entropylow-entropyhigh);
		}
//		std::cout<<"\n";
		ret.p2=this->getP(diffs,ret.diff);
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
