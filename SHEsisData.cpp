/*
 * SHEsisData.cpp
 *
 *  Created on: Aug 3, 2014
 *      Author: ionadmin
 */

#include "SHEsisData.h"

namespace SHEsis {

SHEsisData::SHEsisData(int CaseNum, int ControlNum,  int SnpNum, int NumOfChrSet):
		CaseData(CaseNum,SnpNum,NumOfChrSet),
		ControlData(ControlNum,SnpNum,NumOfChrSet){
	// TODO Auto-generated constructor stub


}

SHEsisData::~SHEsisData() {
	// TODO Auto-generated destructor stub
}

}
