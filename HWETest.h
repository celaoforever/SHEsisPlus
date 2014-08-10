/*
 * HWETest.h
 *
 *  Created on: Aug 10, 2014
 *      Author: ionadmin
 */

#ifndef HWETEST_H_
#define HWETEST_H_
#include "SHEsisData.h"
namespace SHEsis {
struct HWETestResult{
	double PearsonP;
	double FisherP;
	double ChiSquare;
};

class HWETest {
public:
	SHEsisData& data;
	std::vector<HWETestResult> vHWETestResult;
	HWETest(SHEsisData& data);
	virtual ~HWETest();
	void AllSnpHWETest();
private:
	void SingleSnpHWETest(int iSnp, double& chi, double& PearsonP, double& FisherP);
};

} /* namespace SHEsis */

#endif /* HWETEST_H_ */
