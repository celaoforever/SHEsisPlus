/*
 * HWETest.cpp
 *
 *  Created on: Aug 10, 2014
 *      Author: ionadmin
 */
#include <boost/assert.hpp>
#include "HWETest.h"
#include "Multinominal.h"
#include "fisher.h"
#include "utility.h"
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

namespace SHEsis {

HWETest::HWETest(SHEsisData& data): data(data),
		vHWETestResult(data.getSnpNum()){
}

HWETest::~HWETest() {
	this->vHWETestResult.clear();
}

void ResetMapValue(boost::unordered_map<short,size_t>& map){
	boost::unordered_map<short, double> ::iterator map_it;
	for(map_it=map.begin();	map_it != map.end();map_it++){
		map[map_it->first]=0;
	};
}

std::string getStrFromSortedShort(std::vector<size_t> v){
	std::stringstream ss;
	for(int i=0;i<v.size()-1;i++)
	{
		ss<<v[i]<<"/";
	}
	ss<<v[v.size()-1];
	return ss.str();
}

double GetExpectedGenotypeFreq(std::string& genotype,
		boost::unordered_map<short,double>& AlleleCount,
		boost::unordered_map<short, size_t>& AlleleType,
		boost::unordered_map<std::string,size_t>& coefficient,
		int SampleNum, int NumOfChrSet){

	ResetMapValue(AlleleType);
	std::vector<std::string> splitted;
	boost::split(splitted, genotype, boost::is_any_of("/"));
	BOOST_ASSERT(splitted.size()>0);

	boost::unordered_map<short,size_t>::iterator iter;
	for(int i=0;i<splitted.size();i++)
	{
		iter=AlleleType.find(std::atoi(splitted[i].c_str()));
		BOOST_ASSERT(AlleleType.end() != iter);
		AlleleType[iter->first]++;
	}

	std::vector<size_t> exponent;
	for(iter=AlleleType.begin();iter != AlleleType.end();iter++){
		exponent.push_back(iter->second);
	};
	std::sort(exponent.begin(),exponent.end());
	std::string AlleleTypeStatus=getStrFromSortedShort(exponent);

	boost::unordered_map<std::string,size_t>::iterator iter2;
	size_t CurCoefficient;
	if(coefficient.end()!=coefficient.find(AlleleTypeStatus)){
		CurCoefficient=iter2->second;
	}else
	{
		CurCoefficient=multi<size_t>(exponent);
		coefficient[AlleleTypeStatus]=CurCoefficient;
	};

	double res=CurCoefficient;
	for(iter=AlleleType.begin();iter != AlleleType.end();iter++){
		short CurAllele=iter->first;
		size_t CurCount=iter->second;
		if( 0 == CurCount)
			continue;
		res*=pow(AlleleCount[CurAllele]/((double)SampleNum*NumOfChrSet),CurCount);
	}

	splitted.clear();
	exponent.clear();

	return res;
}

void HWETest::SingleSnpHWETest(int iSnp, double& chi, double& PearsonP, double& FisherP){
	BOOST_ASSERT(this->data.vLocusInfo[iSnp].CaseAlleleCount.size()
			== this->data.vLocusInfo[iSnp].ControlAlleleCount.size());

	boost::unordered_map<short, size_t> AlleleType;
	boost::unordered_map<short, double> ::iterator map_it;
	for(map_it=this->data.vLocusInfo[iSnp].CaseAlleleCount.begin();
			map_it != this->data.vLocusInfo[iSnp].CaseAlleleCount.end();
			map_it++){
		AlleleType[map_it->first]=0;
	};



}




} /* namespace SHEsis */
