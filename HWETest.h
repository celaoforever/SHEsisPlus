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

	double ControlPearsonP;
	double ControlChiSquare;

	double BothPearsonP;
	double BothChiSquare;
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
	void SingleSnpHWETest(int iSnp, double& CaseChi, double& CasePearsonP,
			double& ControlChi,double& ControlPearsonP,
			double& BothChi, double& BothPearsonP);
	boost::unordered_map<std::string,size_t> vCoefficient;
};

} /* namespace SHEsis */

#endif /* HWETEST_H_ */
