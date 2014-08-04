/*
 * SHEsisData.h
 *
 *  Created on: Aug 3, 2014
 *      Author: ionadmin
 */

#ifndef SHESISDATA_H_
#define SHESISDATA_H_
#include "GenotypeMatrix.h"
#include <boost/unordered_map.hpp>
#include <string.h>
namespace SHEsis {

class SHEsisData {
public:
	SHEsisData(int CaseNum, int ControlNum, int SnpNum, int NumOfChrSet);
	virtual ~SHEsisData();
	void statAlleleFreq();
	void statGenoFreq();
	GenotypeMatrix const  CaseData;
	GenotypeMatrix const  ControlData;
	boost::unordered_map<short const ,std::string const > Num2Allele;
	boost::unordered_map<std::string const ,short const > Allele2Num;
};

}

#endif /* SHESISDATA_H_ */
