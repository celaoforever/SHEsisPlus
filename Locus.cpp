/*
 * Locus.cpp
 *
 *  Created on: Aug 5, 2014
 *      Author: ada
 */

#include "Locus.h"

namespace SHEsis {

Locus::Locus(int mSampleNum, int mNumOfChrSet):data(boost::extents[mSampleNum][mNumOfChrSet]),
		SampleNum(mSampleNum),NumOfChrSet(mNumOfChrSet)
{


}

Locus::~Locus() {
	// TODO Auto-generated destructor stub
}

} /* namespace SHEsis */
