/*
 * HaplotypeDiploid.h
 *
 *  Created on: Sep 9, 2014
 *      Author: ada
 */

#ifndef HAPLOTYPEDIPLOID_H_
#define HAPLOTYPEDIPLOID_H_
#include "HaplotypeBase.h"

namespace SHEsis {

class HaplotypeDiploid: public HaplotypeBase {
public:
	HaplotypeDiploid(boost::shared_ptr<SHEsisData> data);
	HaplotypeDiploid(boost::shared_ptr<SHEsisData> data, int Snp, std::vector<short> mask);
	virtual ~HaplotypeDiploid();
};

} /* namespace SHEsis */

#endif /* HAPLOTYPEDIPLOID_H_ */
