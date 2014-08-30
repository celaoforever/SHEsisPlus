/*
 * Haplotype_test.cpp
 *
 *  Created on: Aug 19, 2014
 *      Author: ada
 */

#include "Haplotype.h"
#include <iostream>
#include <boost/test/minimal.hpp>
#include <boost/assert.hpp>
#include <initializer_list>

void testHp(){
	int sampleNum=5;
	int snpNum=3;
	int chrSetNum=4;
	SHEsis::SHEsisData data(sampleNum,snpNum,chrSetNum);
	for(int iSample=0;iSample<sampleNum;iSample++){
		data.vLabel[iSample]=SHEsis::CASE;
	}

	//data.mGenotype=
//	int a[5][3][4]=
//	{{{1,1,4,4},{3,3,2,2},{4,4,2,2}},
//	{{1,1,3,3},{3,3,2,2},{4,4,2,2}},
//	{{1,1,3,3},{3,3,3,3},{4,4,2,2}},
//	{{1,1,3,3},{3,3,2,2},{4,4,2,2}},
//	{{1,1,3,3},{3,3,2,2},{2,2,2,2}}};

//	int a[5][3][4]=
//	{{{0,1,4,4},{3,3,2,2},{4,4,2,2}},
//	{{1,1,3,3},{3,3,2,2},{4,4,2,2}},
//	{{1,1,3,3},{3,3,3,3},{4,4,2,2}},
//	{{1,1,3,3},{3,3,2,2},{4,4,2,2}},
//	{{1,1,3,3},{3,3,2,2},{2,2,2,2}}};

	int a[5][3][4]=
	{{{0,1,4,4},{3,3,2,2},{4,4,4,2}},
	{{1,1,1,1},{3,3,2,2},{4,2,2,2}},
	{{1,1,3,3},{3,3,3,3},{4,4,2,2}},
	{{1,1,3,3},{3,3,2,2},{4,4,2,2}},
	{{1,1,3,3},{0,0,0,0},{2,2,2,2}}};
	for(int iSample=0;iSample<sampleNum;iSample++){
		BOOST_ASSERT(iSample<data.vLabel.size());
		data.vLabel[iSample]=SHEsis::CASE;
		for(int iSnp=0;iSnp<snpNum;iSnp++){
			for(int iChrset=0;iChrset<chrSetNum;iChrset++){
				data.mGenotype[iSample][iSnp][iChrset]=a[iSample][iSnp][iChrset];
			}
		}
	}
	std::vector<short> mask(3);
	mask[0]=0;
	mask[1]=1;
	mask[2]=1;
	data.statCount(data.vLabel);
	data.printLocusInfo();
	SHEsis::Haplotype hp(data,2,mask);
	//SHEsis::Haplotype hp(data);
	//hp.statOccurence();
	//SHEsis::IndexingVariables var;
	//SHEsis::IndexingVariables var2=
	hp.BuildModel(1);

}

void testAs(){
	boost::shared_ptr<int[]> sizes(new int[3]);
	sizes[0]=2;
	sizes[1]=1;
	sizes[2]=2;

	boost::shared_ptr<int[]> array(new int[2]);
	array[0]=1;
	array[1]=1;

	SHEsis::ArrayStorage AS1(sizes);
	SHEsis::ArrayStorage AS2(sizes,array);
	//AS2.get(array);
	AS2.getArray();
	AS2.getDimension();
	AS2.getSizes();
	//AS2.set(array,1);
	AS2.set(1,1);
}

void testIV()
{
	SHEsis::IndexingVariables iv;
	iv.add("test",SetSharedPtr(2,1,1));
	iv.getEnumeration("test",SetSharedPtr(1,0));
	iv.printHmKey();
}
int
test_main(int,char*[]){

	testHp();
	return boost::exit_success;
}

