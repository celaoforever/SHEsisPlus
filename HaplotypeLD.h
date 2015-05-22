/*
 * HaplotypeLD.h
 *
 *  Created on: May 22, 2015
 *      Author: ionadmin
 */

#ifndef HAPLOTYPELD_H_
#define HAPLOTYPELD_H_
#include "HaplotypeEM.h"
namespace SHEsis {

class HaplotypeLD {
public:
	HaplotypeLD(boost::shared_ptr<SHEsisData> data);
	HaplotypeLD(boost::shared_ptr<SHEsisData> data, int Snp,
	              std::vector<short> mask);
	virtual ~HaplotypeLD();
private:
	boost::shared_ptr<HaplotypeEM> hapEM;
	boost::shared_ptr<SHEsisData> data;
	std::vector<int> SnpIdx;
	double ldT;
	boost::shared_ptr<double[]> ldPattern;
	std::vector<short> mask;
	void getHaplotypeSub(int start, int end);
	void getAdjacentLD();
	std::vector<int>  getLDBlock(double t);
	boost::shared_ptr<HaplotypeEM> hap;

};

} /* namespace SHEsis */

#endif /* HAPLOTYPELD_H_ */
