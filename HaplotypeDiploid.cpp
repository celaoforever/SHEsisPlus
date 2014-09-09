/*
 * HaplotypeDiploid.cpp
 *
 *  Created on: Sep 9, 2014
 *      Author: ada
 */

#include "HaplotypeDiploid.h"

namespace SHEsis {

HaplotypeDiploid::HaplotypeDiploid(boost::shared_ptr<SHEsisData> data):HaplotypeBase(data) {

}

HaplotypeDiploid::HaplotypeDiploid(boost::shared_ptr<SHEsisData>  data, int Snp, std::vector<short> mask):HaplotypeBase(data,mask){

}

HaplotypeDiploid::~HaplotypeDiploid() {

}

} /* namespace SHEsis */
