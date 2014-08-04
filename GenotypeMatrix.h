/*
 * GenotypeMatrix.h
 *
 *  Created on: Aug 3, 2014
 *      Author: ionadmin
 */

#ifndef GENOTYPEMATRIX_H_
#define GENOTYPEMATRIX_H_
#include <boost/multi_array.hpp>
#include <boost/assert.hpp>
namespace SHEsis {

class GenotypeMatrix {
public:
	GenotypeMatrix(int SampleNum, int SnpNum, int NumOfChrSet);
	short GetAlleleBySampleSnpChrSet(int sample, int snp, int chrset);
	void SetAlleleBySampleSnpChrSet(int sample, int snp, int chrset,short val);
	int getNumOfChrSet(){return this->NumOfChrSet;};
	int getSampleNum(){return this->SampleNum;};
	int getSnpNum(){return this->SnpNum;};
	virtual ~GenotypeMatrix();

private:
	 boost::multi_array<  short, 3> genotype;
	 int SampleNum;
	 int SnpNum;
	 int NumOfChrSet;

};

};

#endif /* GENOTYPEMATRIX_H_ */
