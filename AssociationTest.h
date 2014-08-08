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
typedef enum{
	PERSON,
	FISHER
} TestType;

struct LocusAssiciationTestResult{
	double GenoTypePersonP;
	double GenotypeFisherP;
	double GenotypeChiSquare;

	double AllelePersonP;
	double AlleleFisherP;
	double AleleChiSquare;
	double AlleleOddsRatio;
};

class AssociationTest {
public:
	SHEsisData&  data;
	std::vector<LocusAssiciationTestResult> vAssocationTestResult;
	AssociationTest();
	void SingleSnpTest(int iSnp, TestType method, double &p, int &ndf);
	void
	virtual ~AssociationTest();
private:
	int NumOfPermutation;
	TestType method;

};

} /* namespace SHEsis */

#endif /* ASSOCIATIONTEST_H_ */
