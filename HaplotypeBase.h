/*
 * HaplotypeBase.h
 *
 *  Created on: Sep 9, 2014
 *      Author: ada
 */

#ifndef HAPLOTYPEBASE_H_
#define HAPLOTYPEBASE_H_
#include <boost/shared_ptr.hpp>
#include "SHEsisData.h"
namespace SHEsis {


struct HapTestResult{
	std::vector<boost::shared_ptr<short[]> > haplotypes;
	boost::multi_array<short,2> genotypes;
	boost::shared_ptr<int[]> CaseCount;
	boost::shared_ptr<int[]> ControlCount;
	double FisherP;
	double PearsonP;
	double ChiSquare;
	HapTestResult(int sample,int ploidy):
		genotypes(boost::extents[sample][ploidy]),FisherP(1),PearsonP(1),ChiSquare(0)
	{};
	~HapTestResult(){
		haplotypes.clear();
	}
};



class HaplotypeBase {
public:
	boost::shared_ptr<SHEsisData> data;
	HaplotypeBase(boost::shared_ptr<SHEsisData> data,std::vector<short> mask):data(data),mask(mask),Results(data->getSampleNum(),data->getNumOfChrSet()){};
	HaplotypeBase(boost::shared_ptr<SHEsisData> data):data(data),Results(data->getSampleNum(),data->getNumOfChrSet()){};
	virtual ~HaplotypeBase(){
		this->mask.clear();
		this->SnpIdx.clear();
	}
	virtual void startHaplotypeAnalysis()=0;
protected:
	std::vector<short> mask;
	std::vector<int> SnpIdx;
	HapTestResult Results;

};

} /* namespace SHEsis */

#endif /* HAPLOTYPEBASE_H_ */
