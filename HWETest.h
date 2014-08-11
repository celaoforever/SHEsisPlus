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
	double CasePearsonP;
	double CaseFisherP;
	double CaseChiSquare;

	double ControlPearsonP;
	double ControlFisherP;
	double ControlChiSquare;

	double BothPearsonP;
	double BothFisherP;
	double BothChiSquare;
};
typedef enum{
	CASE_,
	CONTROL_,
	BOTH_
} category;

class HWETest {
public:
	SHEsisData& data;
	std::vector<HWETestResult> vHWETestResult;
	HWETest(SHEsisData& data);
	virtual ~HWETest();
	void AllSnpHWETest();
	void printHWETestResults();
private:
	void SingleSnpHWETest(int iSnp, double& CaseChi, double& CasePearsonP, double& CaseFisherP,
			double& ControlChi,double& ControlPearsonP, double& ControlFisherP,
			double& BothChi, double& BothPearsonP, double& BothFisherP);
	boost::unordered_map<std::string,size_t> vCoefficient;
};

} /* namespace SHEsis */

#endif /* HWETEST_H_ */
