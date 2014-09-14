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
	virtual void startHaplotypeAnalysis(){};
	void initSubArray();
private:
	boost::multi_array<int, 3> PhasedData;
	boost::multi_array<int, 3> CurGenotype; //2 loci, diploid
	std::vector<int> CurGenotypeCount;
	boost::shared_ptr<int> Sample2Genotype;
	int phased;

};

} /* namespace SHEsis */

#endif /* HAPLOTYPEDIPLOID_H_ */
