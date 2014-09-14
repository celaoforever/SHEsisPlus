/*
 * HaplotypeDiploid.cpp
 *
 *  Created on: Sep 9, 2014
 *      Author: ada
 */

#include "HaplotypeDiploid.h"

namespace SHEsis {

HaplotypeDiploid::HaplotypeDiploid(boost::shared_ptr<SHEsisData> data):HaplotypeBase(data),phased(0),CurGenotype(boost::extents[2][2][2]) {
BOOST_ASSERT(data->getNumOfChrSet() == 2);
}

HaplotypeDiploid::HaplotypeDiploid(boost::shared_ptr<SHEsisData>  data, int Snp, std::vector<short> mask):HaplotypeBase(data,mask),phased(0),CurGenotype(boost::extents[2][2][2]){
BOOST_ASSERT(data->getNumOfChrSet() == 2);
}

HaplotypeDiploid::~HaplotypeDiploid() {
}
template <typename T>
bool isMissing(boost::multi_array<T,3> data, int sample, int start, int end){
	BOOST_ASSERT(sample<data.shape()[0]);
	BOOST_ASSERT(start<=end);
	BOOST_ASSERT(end<data.shape()[1]);
	for(int i=start;i<=end;i++){
		for(int j=0;j<data.shape()[2];j++){
			if(0 == data[sample][i][j])
				return true;
		}
	}
	return false;
}
template bool isMissing(boost::multi_array<int,3> data, int sample, int start, int end);
template bool isMissing(boost::multi_array<short,3> data, int sample, int start, int end);
template <typename T>
bool compareGenotype(boost::multi_array<T,3> data, int sample1, int sample2, int start, int end){
	BOOST_ASSERT(sample1<data.shape()[0]);
	BOOST_ASSERT(sample2<data.shape()[0]);
	BOOST_ASSERT(start<=end);
	BOOST_ASSERT(end<data.shape()[1]);
	//for(int )

}

void HaplotypeDiploid::initSubArray(){
	if(0 == this->phased){

	}
}




} /* namespace SHEsis */
