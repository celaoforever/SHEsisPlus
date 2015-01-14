/*
 * GeneInteractionQTL.cpp
 *
 *  Created on: Jan 10, 2015
 *      Author: ionadmin
 */

#include "GeneInteractionQTL.h"

namespace SHEsis {

GeneInteractionQTL::GeneInteractionQTL(boost::shared_ptr<SHEsisData> data):GeneInteraction(data),
		NumBin(10),MaxBin(100),MinBin(3),MinSamplesPerBin(50){
	this->lr.reset(new linear);

}

GeneInteractionQTL::~GeneInteractionQTL() {
	// TODO Auto-generated destructor stub
}

bool sortSampleByQtl(const qtl2sampleIdx& v1, const qtl2sampleIdx& v2){
	return (v1.qtl<v2.qtl);
}
int GeneInteractionQTL::UpdateBins2(std::vector<int>& snp){
	this->bins.clear();
	std::vector<qtl2sampleIdx> validSamples;
	this->GetNonmissingSamples2(snp,validSamples);
	std::sort(validSamples.begin(),validSamples.end(),sortSampleByQtl);
	int SamplesPerBin=validSamples.size()/this->NumBin;
	//adjust according to min samples per bin
	int maxbin=validSamples.size()/this->MinSamplesPerBin;
	if(this->NumBin>maxbin){
		this->NumBin=maxbin;
		SamplesPerBin=this->MinSamplesPerBin;
	}else{
		int rest=validSamples.size()%this->NumBin;
		if(rest>this->MinSamplesPerBin)
			this->NumBin++;
	}
	std::cout<<"NumBin="<<this->NumBin<<",samplesPerbin="<<SamplesPerBin<<"\n";
	if(this->NumBin<this->MinBin)
		return 0;

	for(int i=0;i<this->NumBin-1;i++){
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
	b.start=validSamples[(this->NumBin-1)*SamplesPerBin].qtl;
	b.end=validSamples[validSamples.size()-1].qtl;
	int c=0;
	for(int j=(this->NumBin-1)*SamplesPerBin;j<validSamples.size();j++){
		b.meanqtl+=validSamples[j].qtl;
		b.SampleIdx.push_back(validSamples[j].idx);
		c++;
	}
	b.meanqtl/=(double)c;
	this->bins.push_back(b);
//	int binIdx=0;
//	for(std::list<bin>::iterator iter=this->bins.begin();iter!=this->bins.end();iter++){
//		std::cout<<"bin "<<binIdx<<",start="<<iter->start<<",end="<<iter->end<<",mean="<<iter->meanqtl
//				<<",samples:";
//		for(int i=0;i<iter->SampleIdx.size();i++){
//			std::cout<<iter->SampleIdx[i]<<",";
//		}
//		std::cout<<"\n";
//		binIdx++;
//	}
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
	std::cout<<"snp_set\tnonmissing\tnum_bin\tcoef\tse\tp\tdetail\n";
	for(int i=0;i<this->res.size();i++){
		std::cout<<res[i].snpset<<"\t"<<res[i].nonmissing<<"\t"<<res[i].numBin<<"\t"<<res[i].coef<<"\t"<<res[i].se<<"\t"<<res[i].p<<"\t";
//		for(int j=0;j<res[i].entropy.size();j++){
//			std::cout<<res[i].entropy[j]<<",";
//		}
		std::cout<<"\n";
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
}
double GeneInteractionQTL::getIntervalQTL(bin& b){
//	return (b.start+b.end)/2.;
	return b.meanqtl;
}

void normalizeVector(std::vector<double>& v, std::vector<double>& out){
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
		double val=(v[i]-mean)/stddev;
		out.push_back(val);
	}
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
	ret.nonmissing=this->UpdateBins2(Snp);
	if(ret.nonmissing == 0 ){
		std::cout<<"number of bins is less than min bin, skip...\n";
		return ret;
	};
	//calculate information iteraction within each bins
	for(std::list<bin>::iterator iter=this->bins.begin();iter != this->bins.end();iter++){
		ret.entropy.push_back(this->getInformationInteraction(iter->SampleIdx,Snp));
		ret.qtl.push_back(this->getIntervalQTL(*iter));
	}
	//linear regression
	std::vector<double> norm_qtl;
	std::vector<double> norm_entropy;
	normalizeVector(ret.qtl,norm_qtl);
	normalizeVector(ret.entropy,norm_entropy);
//	lr->resetResponse(ret.qtl);
//	lr->resetSNP(ret.entropy);
	lr->resetResponse(norm_qtl);
	lr->resetSNP(norm_entropy);
	lr->regress();
	for(int i=0;i<ret.qtl.size();i++){
		std::cout<<lr->regressors(0,i)<<"\t\t"<<lr->responses(i)<<"\t\t"<<lr->predictions(i)<<"\n";
	}
	std::cout<<"regression function:\ny="<<lr->coef(1)<<"x"<<(lr->coef(0)>0?"+":"-")<<lr->coef(0)<<",se="<<lr->se(1)<<"\n";
	std::cout<<"se="<<lr->se(1)<<",p="<<lr->p(1)<<"\n\n";

	ret.coef=lr->coef(1);
	ret.p=lr->p(1);
	ret.se=lr->se(1);
	ret.numBin=this->bins.size();
	return ret;
}




} /* namespace SHEsis */
