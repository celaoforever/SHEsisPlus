/*
 * Locus.h
 *
 *  Created on: Aug 5, 2014
 *      Author: ada
 */

#ifndef LOCUS_H_
#define LOCUS_H_
#include <boost/multi_array.hpp>
#include <boost/unordered_map.hpp>
#include <boost/container/vector.hpp>
#include <boost/assert.hpp>
#include <string>

namespace SHEsis {

class Locus {
public:
	Locus(int mSampleNum, int mNumOfChrSet);
	virtual ~Locus();
	double GetFreqByAllele(short a);
	double GetFreqByGenotype(std::string);
	boost::container::vector<double> GetAlleleFreqByGenotype(std::string);
	boost::multi_array<short, 2 > data; //sorted, eg. 0/1/2
	void CalAlleleFreq();
	void CalGenotypeFreq();
private:
	boost::unordered_map<short, std::string> Short2Allele;
	boost::unordered_map< std::string, double> GenotypeFreq; //string is like 1/2/2 coding, sorted
	boost::unordered_map<short, double> AlleleFreq;
	const int SampleNum;
	const int NumOfChrSet;
	std::string LocusName;
};

} /* namespace SHEsis */

#endif /* LOCUS_H_ */
