/*
 * LDTest.cpp
 *
 *  Created on: Sep 17, 2014
 *      Author: ada
 */

#include "LDTest.h"
#include "HaplotypeDiploid.h"
#include "Haplotype.h"
namespace SHEsis {

LDTest::LDTest(boost::shared_ptr<SHEsisData> data):
		data(data),
		res(boost::extents[this->data->getSnpNum()][this->data->getSnpNum()]) {
	// TODO Auto-generated constructor stub

}

LDTest::~LDTest() {
	// TODO Auto-generated destructor stub
}

double LDTest::TwoLociLDTest(int snp1,int snp2){
	double D=0;
	std::vector<int> mask;
	for(int i=0;i<this->data->getSnpNum();i++){
		if(i == snp1 || i == snp2)
			mask.push_back(1);
		else
			mask.push_back(0);
	}
	BOOST_ASSERT(mask.size()==this->data->getSnpNum());
	if(this->data->getNumOfChrSet()<= 2){
		this->hp.reset(new HaplotypeDiploid(this->data,2,mask));
	}else{
		this->hp.reset(new Haplotype(this->data,2,mask));
	}
	hp->startHaplotypeAnalysis();

	return D;
}


void LDTest::AllLociLDtest(){
	for(int i=0;i<this->data->getSnpNum();i++){
		for(int j=i+1;j<this->data->getSnpNum();j++){
			this->res[i][j]=this->TwoLociLDTest(i,j);
		}
	}
}

} /* namespace SHEsis */
