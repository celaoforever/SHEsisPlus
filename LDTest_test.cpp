/*
 * LDTest_test.cpp
 *
 *  Created on: Sep 17, 2014
 *      Author: ada
 */

#include "LDTest.h"

int main(){
	int sampleNum=10;
	int snpNum=5;
	int chrSetNum=2;
	boost::shared_ptr<SHEsis::SHEsisData> data(new SHEsis::SHEsisData(sampleNum,snpNum,chrSetNum));
	SHEsis::LDTest ld(data);
	ld.DrawLDMap();
}
