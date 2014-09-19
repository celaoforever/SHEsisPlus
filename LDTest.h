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
#include "BMP.h"
namespace SHEsis {


typedef enum{
	LD_IN_CASE,
	LD_IN_CTRL,
	LD_IN_BOTH
}LD_TYPE;

class LDTest {
public:
	LDTest(boost::shared_ptr<SHEsisData> data);
	virtual ~LDTest();
	void AllLociLDtest();
	void DrawLDMap();
private:
	double TwoLociLDTest(int snp1,int snp2,LD_TYPE type);
	boost::shared_ptr<SHEsisData> data;
	boost::shared_ptr<HaplotypeBase> hp;
	boost::multi_array<double,2> res;
	BMP* ldmap;
	LD_TYPE ldtype;
};

} /* namespace SHEsis */

#endif /* LDTEST_H_ */