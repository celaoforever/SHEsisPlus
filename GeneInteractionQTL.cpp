/*
 * GeneInteractionQTL.cpp
 *
 *  Created on: Jan 10, 2015
 *      Author: ionadmin
 */

#include "GeneInteractionQTL.h"

namespace SHEsis {

GeneInteractionQTL::GeneInteractionQTL(boost::shared_ptr<SHEsisData> data):GeneInteraction(data),
		MaxBin(100),MinBin(3),MinSamplesPerBin(10){
	this->lr.reset(new linear);

}

GeneInteractionQTL::~GeneInteractionQTL() {
	// TODO Auto-generated destructor stub
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

void GeneInteractionQTL::statIntervalSampleNum(bin& b){
	//[start,end)
	for(int i=0;i<this->data->getSampleNum();i++){
		if(this->data->vQuantitativeTrait[i]>=b.start && this->data->vQuantitativeTrait[i]<b.end){
			b.SampleIdx.push_back(i);
		}
	}
}

bool GeneInteractionQTL::ForwardMergeBins(std::list<bin>::iterator& iter){
		//node pointed by iter will be deleted, and merge into iter+1
		if(iter-1 == this->bins.end())
			return false;
		std::list<bin>::iterator next=iter+1;
		next->start=iter->start;
		next->SampleIdx.insert(next->SampleIdx.begin(),iter->SampleIdx.begin(),iter->SampleIdx.end());
		this->bins.erase(iter);
		iter=next;
	return true;
}

bool GeneInteractionQTL::BackwardMergeBins(std::list<bin>::reverse_iterator& iter){
	//node pointed by iter will be deleted, and merge into iter-1
	if(iter-1 == this->bins.rend())
		return false;
	std::list<bin>::reverse_iterator previous=iter+1;
	previous->end=iter->end;
	previous->SampleIdx.insert(previous->SampleIdx.begin(),iter->SampleIdx.begin(),iter->SampleIdx.end());
	this->bins.erase(--(iter.base()));
	iter=previous;
	return true;
}

void GeneInteractionQTL::print(){
	std::cout<<"snp_set\tnonmissing\tnum_bin\tcoef\tse\tp\n";
	for(int i=0;i<this->res.size();i++){
		std::cout<<res[i].snpset<<"\t"<<res[i].nonmissing<<"\t"<<res[i].numBin<<"\t"<<res[i].coef<<"\t"<<res[i].se<<"\t"<<res[i].p<<"\n";
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
	return (b.start+b.end)/2.;
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
	ret.nonmissing=this->UpdateBins(Snp);
	if(this->bins.size()<this->MinBin){
		std::cout<<"number of bins is less than min bin, skip...\n";
		return ret;
	};
	//calculate information iteraction within each bins
	for(std::list<bin>::iterator iter=this->bins.begin();iter != this->bins.end();iter++){
		ret.entropy.push_back(this->getInformationInteraction(iter->SampleIdx,Snp));
		ret.qtl.push_back(this->getIntervalQTL(*iter));
	}
	//linear regression
	lr->resetResponse(ret.qtl);
	lr->resetSNP(ret.entropy);
	lr->regress();
	ret.coef=lr->coef(1);
	ret.p=lr->p(1);
	ret.se=lr->se(1);
	ret.numBin=this->bins.size();
	return ret;
}




} /* namespace SHEsis */
