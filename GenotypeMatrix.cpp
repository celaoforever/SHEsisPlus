/*
 * GenotypeMatrix.cpp
 *
 *  Created on: Aug 3, 2014
 *      Author: ionadmin
 */

#include "GenotypeMatrix.h"

namespace SHEsis {

GenotypeMatrix::GenotypeMatrix(int SampleNum, int SnpNum, int NumOfChrSet)
:SampleNum(SampleNum),SnpNum(SnpNum),NumOfChrSet(NumOfChrSet),genotype(boost::extents[SampleNum][SnpNum][NumOfChrSet])
 {

};

void GenotypeMatrix::SetAlleleBySampleSnpChrSet(int sample, int snp, int chrset,short val){
	BOOST_ASSERT_MSG(chrset<this->NumOfChrSet,"request NumOfChrSet bigger than max");
	BOOST_ASSERT_MSG(sample<this->SampleNum,"request SampleNum  bigger than max");
	BOOST_ASSERT_MSG(snp<this->SnpNum,"request SnpNum set bigger than max");
	this->genotype[sample][snp][chrset]=val;
};


short GenotypeMatrix::GetAlleleBySampleSnpChrSet(int sample, int snp, int chrset){
	BOOST_ASSERT_MSG(chrset<this->NumOfChrSet,"request NumOfChrSet bigger than max");
	BOOST_ASSERT_MSG(sample<this->SampleNum,"request SampleNum  bigger than max");
	BOOST_ASSERT_MSG(snp<this->SnpNum,"request SnpNum set bigger than max");
	return this->genotype[sample][snp][chrset];
};

GenotypeMatrix::~GenotypeMatrix() {
	// TODO Auto-generated destructor stub
}
;
}
