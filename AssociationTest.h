/*
 * AssociationTest.h
 *
 *  Created on: Aug 8, 2014
 *      Author: ada
 */

#ifndef ASSOCIATIONTEST_H_
#define ASSOCIATIONTEST_H_
#include "SHEsisData.h"
namespace SHEsis {

struct LocusAssiciationTestResult{
	double GenoTypePearsonP;
	double GenotypeFisherP;
	double GenotypeChiSquare;

	double AllelePearsonP;
	double AlleleFisherP;
	double AlleleChiSquare;
	double AlleleOddsRatio;
	double AlleleOddsRatioLowLimit;
	double AlleleOddsRatioUpLimit;
};

class AssociationTest {
public:
	SHEsisData&  data;
	std::vector<LocusAssiciationTestResult> vAssocationTestResult;
	AssociationTest(SHEsisData& data);
	void AssociationTestForAllSnpsAllele();
	void AssociationTestForAllSnpsGenotype();
	virtual ~AssociationTest();
	void printAssociationTestResults();
private:
	int NumOfPermutation;
	void SingleSnpTestAllele(int iSnp, double& FisherP, double& PearsonP,
			double& ChiSquare,double& oddsRatio,double& ORLowLimit, double& ORUpLimit);
	void SingleSnpTestGenotype(int iSnp, double &FisherP, double& PearsonP, double & ChiSquare);
};

} /* namespace SHEsis */

#endif /* ASSOCIATIONTEST_H_ */
