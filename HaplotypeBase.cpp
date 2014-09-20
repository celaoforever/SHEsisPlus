/*
 * HaplotypeBase.cpp
 *
 *  Created on: Sep 20, 2014
 *      Author: ionadmin
 */
#include "HaplotypeBase.h"
namespace SHEsis{

void HaplotypeBase::AssociationTest(){

	int haploNum=this->Results.haplotypes.size();
//    std::cout<<"contigency:\n";
    double* contigency= new double[2*haploNum];
    int idx=0;
    for(int i=0;i<haploNum;i++){
    	contigency[idx++]=this->Results.ControlCount[i];
    	contigency[idx++]=this->Results.CaseCount[i];
//    	std::cout<<this->Results.ControlCount[i]<<","<<this->Results.CaseCount[i]<<"\n";
    };

     int nrow=2;
	  double expect = -1.0;
	  double percnt = 100.0;
	  double emin = 0;
	  double pre = 0, prt = 0;
	  int ws = 300000;
	  try{
		  fexact(&nrow, &haploNum, contigency, &nrow, &expect, &percnt, &emin, &prt, &pre, &ws);
		  this->Results.FisherP=pre;
	  }catch(std::runtime_error &){
		  this->Results.FisherP=-1;
	  }
	  //Pearson's ChiSquare test
	  PearsonChiSquareTest(contigency,nrow,haploNum,this->Results.ChiSquare,this->Results.PearsonP);
	  delete[] contigency;
//	  std::cout<<"fisherp:"<<this->Results.FisherP;
//	  std::cout<<"\npearsonp:"<<this->Results.PearsonP;
}
}


