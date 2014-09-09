/*
 * HWETest.h
 *
 *  Created on: Aug 10, 2014
 *      Author: ionadmin
 */

#ifndef HWETEST_H_
#define HWETEST_H_
#include "SHEsisData.h"
#include <boost/shared_ptr.hpp>
namespace SHEsis {
struct HWETestResult{
	double CasePearsonP;
	double CaseChiSquare;
	double CaseFisherP;
	double ControlPearsonP;
	double ControlChiSquare;
	double ControlFisherP;
	double BothPearsonP;
	double BothChiSquare;
	double BothFisherP;
};
typedef enum{
	CASE_,
	CONTROL_,
	BOTH_
} category;

class HWETest {
public:
	boost::shared_ptr<SHEsisData> data;
	std::vector<HWETestResult> vHWETestResult;
	HWETest(boost::shared_ptr<SHEsisData>  data);
	virtual ~HWETest();
	void AllSnpHWETest();
	void printHWETestResults();
private:
	void SingleSnpHWETest(int iSnp, double& CaseChi, double& CasePearsonP, double& CaseFisherP,
			double& ControlChi,double& ControlPearsonP, double & ControlFisherP,
			double& BothChi, double& BothPearsonP, double& BothFisherP);
	boost::unordered_map<std::string,size_t> vCoefficient;
};

} /* namespace SHEsis */

#endif /* HWETEST_H_ */
