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
	this->vCoefficient.clear();
}

void ResetMapValue(boost::unordered_map<short,size_t>& map){
	boost::unordered_map<short, size_t> ::iterator map_it;
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

double GetExpectedGenotypeFreq(std::string genotype,
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

	boost::unordered_map<std::string,size_t>::iterator iter2=coefficient.find(AlleleTypeStatus);
	size_t CurCoefficient;
	if(coefficient.end()!=iter2){
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

void HWETest::SingleSnpHWETest(int iSnp, double& CaseChi, double& CasePearsonP,
		double& ControlChi,double& ControlPearsonP,
		double& BothChi, double& BothPearsonP){
	BOOST_ASSERT(this->data.vLocusInfo[iSnp].CaseAlleleCount.size()
			== this->data.vLocusInfo[iSnp].ControlAlleleCount.size());

	boost::unordered_map<short, size_t> AlleleType;
	boost::unordered_map<short, double> ::iterator map_it;
	for(map_it=this->data.vLocusInfo[iSnp].CaseAlleleCount.begin();
			map_it != this->data.vLocusInfo[iSnp].CaseAlleleCount.end();
			map_it++){
		AlleleType[map_it->first]=0;
	};

	//calculate HWE for cases
	boost::unordered_map<std::string,double>::iterator genotype_iter;
	int NumOfRow=2;//1st row is for observed freq, 2nd row is for expected frequency
	int NumOfCol=this->data.vLocusInfo[iSnp].CaseGenotypeCount.size();
	int totalGenotype=0;
	double* contigency=new double[NumOfRow*NumOfCol];
	int idx=0;
	for(genotype_iter=this->data.vLocusInfo[iSnp].CaseGenotypeCount.begin();
			genotype_iter != this->data.vLocusInfo[iSnp].CaseGenotypeCount.end();
			genotype_iter++){
			double expectedFreq=GetExpectedGenotypeFreq(genotype_iter->first,
					this->data.vLocusInfo[iSnp].CaseAlleleCount,
					AlleleType,
					this->vCoefficient,
					this->data.getCaseNum(),
					this->data.getNumOfChrSet());
			BOOST_ASSERT(idx<NumOfRow*NumOfCol);
			contigency[idx++]=genotype_iter->second;//(this->data.getCaseNum()*this->data.getNumOfChrSet());
			BOOST_ASSERT(idx<NumOfRow*NumOfCol);
			contigency[idx++]=expectedFreq;
			totalGenotype+=genotype_iter->second;
	};

	//Fisher's exact test:
	  double expect = -1.0;
	  double percnt = 100.0;
	  double emin = 0;
	  double pre = 0, prt = 0;
	  int ws = 300000;
	  CaseChi=0;
	  for(int i=0;i<NumOfCol*NumOfRow;i=i+2){
		  if(contigency[i]<1)
			  contigency[i]*=totalGenotype;
		  if(contigency[i+1]<1)
			  contigency[i+1]*=totalGenotype;
		  CaseChi+=(contigency[i+1]-contigency[i])
				  *(contigency[i+1]-contigency[i])/contigency[i+1];
	  }
	  boost::math::chi_squared dist(1);
	  CasePearsonP= boost::math::cdf(boost::math::complement(dist,CaseChi));
	  //fexact(&NumOfRow, &NumOfCol, contigency, &NumOfRow, &expect, &percnt, &emin, &prt, &pre, &ws);
	//  CaseFisherP=pre;
	  //Pearson's ChiSquare test
	  //PearsonChiSquareTest(contigency,NumOfRow,NumOfCol,CaseChi,CasePearsonP);
	  delete[] contigency;
	  contigency = 0 ;

	//calculate HWE for control
	totalGenotype=0;
	NumOfRow=2;//1st row is for observed freq, 2nd row is for expected frequency
	NumOfCol=this->data.vLocusInfo[iSnp].ControlGenotypeCount.size();
	contigency=new double[NumOfRow*NumOfCol];
	idx=0;
	for(genotype_iter=this->data.vLocusInfo[iSnp].ControlGenotypeCount.begin();
			genotype_iter != this->data.vLocusInfo[iSnp].ControlGenotypeCount.end();
			genotype_iter++){
			double expectedFreq=GetExpectedGenotypeFreq(genotype_iter->first,
					this->data.vLocusInfo[iSnp].ControlAlleleCount,
					AlleleType,
					this->vCoefficient,
					this->data.getControlNum(),
					this->data.getNumOfChrSet());
			BOOST_ASSERT(idx<NumOfRow*NumOfCol);
			contigency[idx++]=genotype_iter->second;//(this->data.getControlNum()*this->data.getNumOfChrSet());
			BOOST_ASSERT(idx<NumOfRow*NumOfCol);
			contigency[idx++]=expectedFreq;
			totalGenotype+=genotype_iter->second;
	};

	//Fisher's exact test:
	ControlChi=0;
	for(int i=0;i<NumOfCol*NumOfRow;i=i+2){
		if(contigency[i]<1)
			contigency[i]*=totalGenotype;
		if(contigency[i+1]<1)
			contigency[i+1]*=totalGenotype;
		  ControlChi+=(contigency[i+1]-contigency[i])
				  *(contigency[i+1]-contigency[i])/contigency[i+1];
	}
	ControlPearsonP= boost::math::cdf(boost::math::complement(dist,ControlChi));
	//PearsonChiSquareTest(contigency,NumOfRow,NumOfCol,ControlChi,ControlPearsonP);
	//fexact(&NumOfRow, &NumOfCol, contigency, &NumOfRow, &expect, &percnt, &emin, &prt, &pre, &ws);
	//ControlFisherP=pre;

	//Pearson's ChiSquare test

	delete[] contigency;
	contigency = 0 ;

	//calculate HWE for case and control
	totalGenotype=0;
	NumOfRow=2;//1st row is for observed freq, 2nd row is for expected frequency
	NumOfCol=this->data.vLocusInfo[iSnp].BothGenotypeCount.size();
	contigency=new double[NumOfRow*NumOfCol];
	idx=0;
	for(genotype_iter=this->data.vLocusInfo[iSnp].BothGenotypeCount.begin();
			genotype_iter != this->data.vLocusInfo[iSnp].BothGenotypeCount.end();
			genotype_iter++){
			double expectedFreq=GetExpectedGenotypeFreq(genotype_iter->first,
					this->data.vLocusInfo[iSnp].BothAlleleCount,
					AlleleType,
					this->vCoefficient,
					this->data.getSampleNum(),
					this->data.getNumOfChrSet());
			BOOST_ASSERT(idx<NumOfRow*NumOfCol);
			contigency[idx++]=genotype_iter->second;//(this->data.getSampleNum()*this->data.getNumOfChrSet());
			BOOST_ASSERT(idx<NumOfRow*NumOfCol);
			contigency[idx++]=expectedFreq;
			totalGenotype+=genotype_iter->second;
	};

	//Fisher's exact test:
	BothChi=0;
	for(int i=0;i<NumOfCol*NumOfRow;i=i+2){
		if(contigency[i]<1)
			contigency[i]*=totalGenotype;
		if(contigency[i+1]<1)
			contigency[i+1]*=totalGenotype;
		  BothChi+=(contigency[i+1]-contigency[i])
				  *(contigency[i+1]-contigency[i])/contigency[i+1];
	}
	BothPearsonP= boost::math::cdf(boost::math::complement(dist,BothChi));
	//fexact(&NumOfRow, &NumOfCol, contigency, &NumOfRow, &expect, &percnt, &emin, &prt, &pre, &ws);
	//BothFisherP=pre;
	//Pearson's ChiSquare test
	//PearsonChiSquareTest(contigency,NumOfRow,NumOfCol,BothChi,BothPearsonP);
	delete[] contigency;
	contigency = 0 ;
}

void HWETest::AllSnpHWETest(){
	for(int i=0;i<this->data.getSnpNum();i++){
		this->SingleSnpHWETest(i,this->vHWETestResult[i].CaseChiSquare,this->vHWETestResult[i].CasePearsonP,
				this->vHWETestResult[i].ControlChiSquare,this->vHWETestResult[i].ControlPearsonP,
				this->vHWETestResult[i].BothChiSquare,this->vHWETestResult[i].BothPearsonP);
	};
};

void HWETest::printHWETestResults(){
	for(int i=0;i<this->vHWETestResult.size();i++){
		std::cout<<"\nLocus "<<i<<":\nHWE test for case:\n(chi,pearsonp)=("<<
				this->vHWETestResult[i].CaseChiSquare<<","<<
				this->vHWETestResult[i].CasePearsonP<<")\n"<<
				"HWE test for control:\n(chi,pearsonp)=("<<
				this->vHWETestResult[i].ControlChiSquare<<","<<
				this->vHWETestResult[i].ControlPearsonP<<")\n"<<
				"HWE test for case and control:\n(chi,pearsonp)=("<<
				this->vHWETestResult[i].BothChiSquare<<","<<
				this->vHWETestResult[i].BothPearsonP<<")\n";

	}
}

} /* namespace SHEsis */
