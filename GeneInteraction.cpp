/*
 * GeneInteraction.cpp
 *
 *  Created on: Dec 30, 2014
 *      Author: ionadmin
 */

#include "GeneInteraction.h"
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <boost/assert.hpp>
namespace SHEsis {


genotypeRetCode GeneInteraction::genotypeEqual(std::string geno,int sample,int snp){
	 if(this->data->mGenotype[sample][snp][0]==0)
		 return MISS;
	 std::stringstream str("");//=this->mGenotypeStr[sample][snp];
	 str<<this->data->mGenotype[sample][snp][0];
	 for(int p=1;p<this->data->getNumOfChrSet();p++){
		 str<<"/"<<this->data->mGenotype[sample][snp][p];
	 }
	 return (!std::strcmp(str.str().data(), geno.data())?EQUAL:NOT_EQUAL);
}

GeneInteraction::GeneInteraction(boost::shared_ptr<SHEsisData> data):lowbound(2),upperbound(2),data(data),
		mGenotypeStr(boost::extents[data->getSampleNum()][data->getSnpNum()]){
	if(this->data->vLocusInfo[0].BothGenotypeCount.begin() == this->data->vLocusInfo[0].BothGenotypeCount.end()){
		this->data->statCount(this->data->vLabel);
	}
	for(int sample=0;sample<this->data->getSampleNum();sample++){
		for(int snp=0;snp<this->data->getSnpNum();snp++){
			std::stringstream geno("");
			geno<<(this->data->mGenotype[sample][snp][0]);
			for(int p=1;p<this->data->getNumOfChrSet();p++){
				geno<<"/"<<this->data->mGenotype[sample][snp][p];
			}
			this->mGenotypeStr[sample][snp]=geno.str();
		}
	}


}

GeneInteraction::~GeneInteraction() {
	// TODO Auto-generated destructor stub
}

} /* namespace SHEsis */
