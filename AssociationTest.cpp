/*
 * AssociationTest.cpp
 *
 *  Created on: Aug 8, 2014
 *      Author: ada
 */

#include "AssociationTest.h"
#include <boost/assert.hpp>
#include "utility.h"
#include "fisher.h"
#include <iostream>
#include <math.h>
#include <algorithm>
namespace SHEsis {

AssociationTest::AssociationTest(boost::shared_ptr<SHEsisData>  mdata):data(mdata),
		vAssocationTestResult(mdata->getSnpNum()),
		NumOfPermutation(1000)
{
}

AssociationTest::~AssociationTest() {
	vAssocationTestResult.clear();
}

void AssociationTest::AssociationTestForAllSnpsAllele(){
	for(int i=0;i<this->vAssocationTestResult.size();i++){
		this->SingleSnpTestAllele(i,
				this->vAssocationTestResult[i].AlleleFisherP,
				this->vAssocationTestResult[i].AllelePearsonP,
				this->vAssocationTestResult[i].AlleleChiSquare,
				this->vAssocationTestResult[i].AlleleOddsRatio,
				this->vAssocationTestResult[i].AlleleOddsRatioLowLimit,
				this->vAssocationTestResult[i].AlleleOddsRatioUpLimit);
	}
}

void AssociationTest::AssociationTestForAllSnpsGenotype(){
	for(int i=0;i<this->vAssocationTestResult.size();i++){
		this->SingleSnpTestGenotype(i,
				this->vAssocationTestResult[i].GenotypeFisherP,
				this->vAssocationTestResult[i].GenoTypePearsonP,
				this->vAssocationTestResult[i].GenotypeChiSquare);
	}
}

void getTheSmallestP(std::vector<LocusAssiciationTestResult> res,double& allelep,double& genop){
	allelep=1;
	genop=1;
	for(int i=0;i<res.size();i++){
		if(res[i].AllelePearsonP < allelep){
			allelep=res[i].AllelePearsonP;
		};
		if(res[i].GenoTypePearsonP< genop){
			genop=res[i].GenoTypePearsonP;
		}
	}
}

int getRank(double p, std::vector<double> v){
	for(int i=0;i<v.size()-1;i++){
		if(p>=v[i] && p<=v[i+1])
			return i;
	}
	return v.size();
}

void AssociationTest::permutation(){
	this->vPermutateLabel=this->data->vLabel;
	for(int i=0;i<this->NumOfPermutation;i++){
		std::random_shuffle(this->vPermutateLabel.begin(), this->vPermutateLabel.end());
		this->data->statCount(this->vPermutateLabel);
		this->AssociationTestForAllSnpsAllele();
		this->AssociationTestForAllSnpsGenotype();
		double ap,gp;
		getTheSmallestP(this->vAssocationTestResult,ap,gp);
		this->PermutationPAllele.push_back(ap);
		this->PermutationPGenotype.push_back(gp);
	};
	std::sort (this->PermutationPAllele.begin(), this->PermutationPAllele.end());
	std::sort (this->PermutationPGenotype.begin(), this->PermutationPGenotype.end());
//	std::cout<<"PermutationPAllele:\n";
//	for(int i=0;i<this->PermutationPAllele.size();i++){
//		std::cout<<this->PermutationPAllele[i]<<",";
//	}
//	std::cout<<"\nPermutationPGenotype:\n";
//	for(int i=0;i<this->PermutationPGenotype.size();i++){
//		std::cout<<this->PermutationPGenotype[i]<<",";
//	}

	this->data->statCount(data->vLabel);
	this->AssociationTestForAllSnpsAllele();
	this->AssociationTestForAllSnpsGenotype();
	for(int i=0;i<this->vAssocationTestResult.size();i++){
		this->vAssocationTestResult[i].AllelePermutationP=(double)getRank(
				this->vAssocationTestResult[i].AllelePearsonP,this->PermutationPAllele)
				/(double)this->NumOfPermutation;
		this->vAssocationTestResult[i].GenotypePermutationP=(double)getRank(
				this->vAssocationTestResult[i].GenoTypePearsonP,this->PermutationPGenotype)
				/(double)this->NumOfPermutation;
	};
}


void AssociationTest::printAssociationTestResults()
{
	for(int i=0;i<this->vAssocationTestResult.size();i++){
		std::cout<<"\nLocus "<<i<<"\nAllele association test:\n"<<
				"(fisher's p, pearson's p, chi suqare, permutation p,or, orLow, orUp)="<<
				"("<<vAssocationTestResult[i].AlleleFisherP<<","<<
				vAssocationTestResult[i].AllelePearsonP<<","<<
				vAssocationTestResult[i].AlleleChiSquare<<","<<
				vAssocationTestResult[i].AllelePermutationP<<","<<
				vAssocationTestResult[i].AlleleOddsRatio<<","<<
				vAssocationTestResult[i].AlleleOddsRatioLowLimit<<","<<
				vAssocationTestResult[i].AlleleOddsRatioUpLimit<<")\n";
		std::cout<<"Genotype association test:\n"<<
				"(fisher's p, pearson's p, chi square,permutation p)=("<<
				vAssocationTestResult[i].GenotypeFisherP<<","<<
				vAssocationTestResult[i].GenoTypePearsonP<<","<<
				vAssocationTestResult[i].GenotypeChiSquare<<","<<
				vAssocationTestResult[i].GenotypePermutationP<<")\n";

	}
}


void AssociationTest::SingleSnpTestAllele(int iSnp, double& FisherP, double& PearsonP,
		double& ChiSquare,double& oddsRatio,double& ORLowLimit, double& ORUpLimit){

	BOOST_ASSERT(this->data->vLocusInfo[iSnp].CaseAlleleCount.size()
			== this->data->vLocusInfo[iSnp].ControlAlleleCount.size());

	int NumOfAlleleType=this->data->vLocusInfo[iSnp].ControlAlleleCount.size();
	int NumOfPhenotype=2;
	double* contigency=new double[NumOfPhenotype*NumOfAlleleType];

	boost::unordered_map<short, double> ::iterator map_it;
	int idx=0;
	for(map_it=this->data->vLocusInfo[iSnp].CaseAlleleCount.begin();
			map_it != this->data->vLocusInfo[iSnp].CaseAlleleCount.end();
			map_it++){
		BOOST_ASSERT(this->data->vLocusInfo[iSnp].ControlAlleleCount.end()!=
				this->data->vLocusInfo[iSnp].ControlAlleleCount.find(map_it->first));
		contigency[idx++]=this->data->vLocusInfo[iSnp].ControlAlleleCount[map_it->first];
		contigency[idx++]=map_it->second;
	};

	//Fisher's exact test:
	  double expect = -1.0;
	  double percnt = 100.0;
	  double emin = 0;
	  double pre = 0, prt = 0;
	  int ws = 300000;
	  fexact(&NumOfPhenotype, &NumOfAlleleType, contigency, &NumOfPhenotype, &expect, &percnt, &emin, &prt, &pre, &ws);
	  FisherP=pre;

	  //Pearson's ChiSquare test
	  PearsonChiSquareTest(contigency,NumOfPhenotype,NumOfAlleleType,ChiSquare,PearsonP);

	  //get odds ratio
	  if((2 == NumOfAlleleType)&&( 0 != contigency[3] &&
			  0!= contigency[0] && 0 !=contigency[1] && 0 != contigency[2])){
		  oddsRatio=(contigency[1]*contigency[2]/(contigency[0]*contigency[3]));
		  double v=1/contigency[1]+1/contigency[2]+1/contigency[3]+1/contigency[0];
		  ORLowLimit=oddsRatio*exp(-1.96*sqrt(v));
		  ORUpLimit=oddsRatio*exp(1.96*sqrt(v));
	  }else
	  {
		  oddsRatio=-1;
		  ORLowLimit=-1;
		  ORUpLimit=-1;
	  }
	  delete[] contigency;
	  contigency = 0;
}


void AssociationTest::SingleSnpTestGenotype(int iSnp, double& FisherP, double& PearsonP, double& ChiSquare){

	BOOST_ASSERT(this->data->vLocusInfo[iSnp].CaseGenotypeCount.size()
			== this->data->vLocusInfo[iSnp].ControlGenotypeCount.size());

	int NumOfAlleleType=this->data->vLocusInfo[iSnp].ControlGenotypeCount.size();
	int NumOfPhenotype=2;
	double* contigency=new double[NumOfPhenotype*NumOfAlleleType];

	boost::unordered_map<std::string, double> ::iterator map_it;
	int idx=0;
	for(map_it=this->data->vLocusInfo[iSnp].CaseGenotypeCount.begin();
			map_it != this->data->vLocusInfo[iSnp].CaseGenotypeCount.end();
			map_it++){
		BOOST_ASSERT(this->data->vLocusInfo[iSnp].ControlGenotypeCount.end()!=
				this->data->vLocusInfo[iSnp].ControlGenotypeCount.find(map_it->first));
		contigency[idx++]=this->data->vLocusInfo[iSnp].ControlGenotypeCount[map_it->first];
		contigency[idx++]=map_it->second;
	};

	//Fisher's exact test:
	  double expect = -1.0;
	  double percnt = 100.0;
	  double emin = 0;
	  double pre = 0, prt = 0;
	  int ws = 300000;
	  try{
	  fexact(&NumOfPhenotype, &NumOfAlleleType, contigency, &NumOfPhenotype, &expect, &percnt, &emin, &prt, &pre, &ws);
	  FisherP=pre;
	  }catch(std::runtime_error &){
		  FisherP=-1;
	  }

	  //Pearson's ChiSquare test
	  PearsonChiSquareTest(contigency,NumOfPhenotype,NumOfAlleleType,ChiSquare,PearsonP);
	  delete[] contigency;
	  contigency = 0 ;
}


} /* namespace SHEsis */
