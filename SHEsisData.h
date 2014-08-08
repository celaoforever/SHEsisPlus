/*
 * SHEsisData.h
 *
 *  Created on: Aug 3, 2014
 *      Author: ionadmin
 */

#ifndef SHESISDATA_H_
#define SHESISDATA_H_
#include <boost/multi_array.hpp>
#include <boost/unordered_map.hpp>
#include <vector>
#include <string>

#define GENOTYPE_MISSING 0
namespace SHEsis {

struct LocusInfo{
	std::string name;
	boost::unordered_map<short, double> CaseAlleleCount;
	boost::unordered_map<short, double> ControlAlleleCount;
	boost::unordered_map<std::string, double> CaseGenotypeCount;
	boost::unordered_map<std::string, double> ControlGenotypeCount;
};

typedef enum{
	MISSING,
	CONTROL,
	CASE
} SampleStatus;

class SHEsisData {
public:
	SHEsisData(int SampleNum, int SnpNum, int NumOfChrSet);
	int getNumOfChrSet(){return this->NumOfChrSet;};
	int getSampleNum(){return this->SampleNum;};
	int getSnpNum(){return this->SnpNum;};
	virtual ~SHEsisData();
	boost::multi_array< short, 3> mGenotype;
	std::vector< SampleStatus > vLabel;
	std::vector<SampleStatus> vPermutateLabel;
	void statCount(std::vector< SampleStatus > & label);
	void printLocusInfo();
private:
	boost::unordered_map<short, std::string> code2allele;
	std::vector<LocusInfo> vLocusInfo;
	const int SampleNum;
	const int SnpNum;
	const int NumOfChrSet;

};

};


#endif /* SHESISDATA_H_ */
