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

struct HaploPair{
	int hap1;
	int hap2;
	bool operator==(const HaploPair& p){
		return(((hap1==p.hap1)&&(hap2==p.hap2))||
				((hap1==p.hap2)&&(hap2==p.hap1)));
	}
};

struct OneGenotypeExpandedHaplo{
	  static std::vector< boost::shared_ptr<short[]> > haploType;
	  static std::vector<double> hapfreq;
	  std::vector<HaploPair> hp;
};


class HaplotypeDiploid: public HaplotypeBase {
public:
	HaplotypeDiploid(boost::shared_ptr<SHEsisData> data);
	HaplotypeDiploid(boost::shared_ptr<SHEsisData> data, int Snp, std::vector<short> mask);
	virtual ~HaplotypeDiploid();
	virtual void startHaplotypeAnalysis(){};
	void GenerateUniqueGenotype();
	void GenerateInterMediate();
	void ReturnGenotypeCode(int sample,short& geno1, short& geno2);
	OneGenotypeExpandedHaplo OneGenoExpandHaplo(int sample);
	void ExpandAllGenotype();
	void CalculateFreq();
	void generateAllPossibleHap();
private:
	boost::multi_array<short, 3> PhasedData;
	std::vector<int> UniqueGenotypeIdx; //2 loci, diploid
	std::vector<int> UniqueGenotypeCount;
	boost::multi_array<short, 3> InterMediate;
	std::vector<std::string> InterMediateGenoCode;
	std::vector<int> CurGenotypeCount;
	boost::shared_ptr<int> Sample2Genotype;
	std::vector<OneGenotypeExpandedHaplo> Expanded;
	int phased;
	double err;

};

} /* namespace SHEsis */

#endif /* HAPLOTYPEDIPLOID_H_ */
