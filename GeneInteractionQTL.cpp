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
	// TODO Auto-generated constructor stub

}

GeneInteractionQTL::~GeneInteractionQTL() {
	// TODO Auto-generated destructor stub
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

gxgQTLRes GeneInteractionQTL::GetOneSNPCombinationInformationGain(std::vector<int>& Snp){

}



} /* namespace SHEsis */
