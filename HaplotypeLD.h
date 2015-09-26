/*
 * HaplotypeLD.h
 *
 *  Created on: May 22, 2015
 *      Author: ionadmin
 */

#ifndef HAPLOTYPELD_H_
#define HAPLOTYPELD_H_
#include "HaplotypeEM.h"
#include <boost/bimap.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/tags/tagged.hpp>
namespace SHEsis {
struct haplo{};
struct code{};
typedef boost::bimap<boost::bimaps::tagged<std::string,haplo>,boost::bimaps::tagged<short,code> > hapMap;
class HaplotypeLD {
public:
	HaplotypeLD(boost::shared_ptr<SHEsisData> data);
	HaplotypeLD(boost::shared_ptr<SHEsisData> data, int Snp,
	              std::vector<short> mask);
	void phaseAll();
	void getResults();
	virtual ~HaplotypeLD();
	std::string reporthtml();
	std::string reporttxt();
	void AssociationTest();
	void setAdjust(bool t){this->adjust=t;};
	void setSilent(bool t){this->silent=t;};
	std::vector<boost::shared_ptr<short[]> > haplotypes;
	void setLDT(double t){this->ldT=t;};
	void setFreqThreshold(double t){this->lft=t;};
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
	void reducePhasedHap(boost::multi_array<short, 3> PhasedData,int curBlockIdx);
	boost::shared_ptr<HaplotypeEM> hap;
	boost::shared_ptr<SHEsisData> blockData;
	hapMap hapcode;
	short curCode;
	bool adjust;
	bool silent;
	double lft;
};

} /* namespace SHEsis */

#endif /* HAPLOTYPELD_H_ */
