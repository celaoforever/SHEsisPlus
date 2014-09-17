/*
 * LDTest.h
 *
 *  Created on: Sep 17, 2014
 *      Author: ada
 */

#ifndef LDTEST_H_
#define LDTEST_H_
#include "HaplotypeBase.h"
#include <boost/shared_ptr.hpp>
namespace SHEsis {

class LDTest {
public:
	LDTest(boost::shared_ptr<SHEsisData> data);
	virtual ~LDTest();
	void AllLociLDtest();
private:
	double TwoLociLDTest(int snp1,int snp2);
	boost::shared_ptr<SHEsisData> data;
	boost::shared_ptr<HaplotypeBase> hp;
	boost::multi_array<double,2> res;
};

} /* namespace SHEsis */

#endif /* LDTEST_H_ */
