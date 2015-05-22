/*
 * HaplotypeLD.cpp
 *
 *  Created on: May 22, 2015
 *      Author: ionadmin
 */

#include "HaplotypeLD.h"
#include "LDTest.h"
namespace SHEsis {

HaplotypeLD::HaplotypeLD(boost::shared_ptr<SHEsisData> data):data(data){
	  for (int i = 0; i < this->data->getSnpNum(); i++) {
	    this->SnpIdx.push_back(i);
	  }
};

HaplotypeLD::HaplotypeLD(boost::shared_ptr<SHEsisData> data, int Snp,
              std::vector<short> mask){
	  for (int i = 0; i < mask.size(); i++) {
	    if (mask[i]) {
	      this->SnpIdx.push_back(i);
	    };
	  };
	  this->mask=mask;
};

void HaplotypeLD::getAdjacentLD(){
	  LDTest ld(this->data);
	  if (this->data->vLocusInfo[0].BothAlleleCount.size() == 0 &&
	      this->data->vQuantitativeTrait.size() > 0) {
	    this->data->statCount();
	  } else if (this->data->vLocusInfo[0].BothAlleleCount.size() == 0 &&
	             this->data->vQuantitativeTrait.size() == 0)
	    this->data->statCount(this->data->vLabel);
	  this->ldPattern.reset(new double[this->SnpIdx.size()-1]);
	  for(int i=0;i<this->SnpIdx.size()-1;i++){
		  double R2,D;
		  ld.TwoLociLDTest(SnpIdx[i], SnpIdx[i+1], LD_IN_BOTH, R2, D);
		  this->ldPattern[i]=R2;
	  }
}

std::vector<int> HaplotypeLD::getLDBlock(double t){
	this->getAdjacentLD();
	std::vector<int> ret;
	ret.push_back(0);
	for(int i=0;i<this->data->getSnpNum()-1;i++){
		if(this->ldPattern[i]<t){
			ret.push_back(i+1);
		}
	}
	ret.push_back(this->SnpIdx.size());
	//snp from ret[i] and ret[i+1]-1 is a block
	return ret;
}


void HaplotypeLD::getHaplotypeSub(int start, int end){
	std::vector<short>::const_iterator first=this->mask.begin()+this->SnpIdx[start];
	std::vector<short>::const_iterator last=this->mask.begin()+this->SnpIdx[end]+1;
	std::vector<short> m(first,last);
	this->hap.reset(new HaplotypeEM(this->data,end-start+1,m));
	this->hap->startHaplotypeAnalysis();
}

HaplotypeLD::~HaplotypeLD() {
	// TODO Auto-generated destructor stub
}

} /* namespace SHEsis */
