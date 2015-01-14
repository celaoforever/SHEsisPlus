/*
 * GeneInteractionQTL_test.cpp
 *
 *  Created on: Jan 10, 2015
 *      Author: ionadmin
 */

#include "GeneInteractionQTL.h"

using namespace std;
using namespace SHEsis;
#include <boost/random/discrete_distribution.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
boost::mt19937 boost_rng;


boost::shared_ptr<SHEsis::SHEsisData> GenerateRandomData(int sampleNum,
                                                         int snpNum,
                                                         int chrSetNum) {
  double GenoProb[] = {0.0, 0.5,0.5};  // 0.01 is missing genotype for individuals
  boost::normal_distribution<> nd(0,1); //phenotype, sigma=2, mean=10
  boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > var_nor(boost_rng, nd);
  boost::random::discrete_distribution<> distGeno(GenoProb);
  boost::shared_ptr<SHEsis::SHEsisData> data(
      new SHEsis::SHEsisData(sampleNum, snpNum, chrSetNum));
  for (int iSample = 0; iSample < sampleNum; iSample++) {
	  data->vQuantitativeTrait.push_back(var_nor());
    for (int iSnp = 0; iSnp < snpNum; iSnp++) {
      for (int iChrset = 0; iChrset < chrSetNum; iChrset++) {
        data->mGenotype[iSample][iSnp][iChrset] = distGeno(boost_rng);
      }
    }
  }
  return data;
}

void SetGXG(boost::shared_ptr<SHEsis::SHEsisData> data,int snp1,int snp2,double pre){
	int ploidy=data->getNumOfChrSet();
	boost::uniform_int<> dist(1,100);
	boost::variate_generator<boost::mt19937,boost::uniform_int<> > dice(boost_rng,dist);
	for(int i=0;i<data->getSampleNum();i++){
		double n=(double)dice()/100.;
		if(n>pre){
			continue;
		}
		double qtl=data->vQuantitativeTrait[i];
		for(int p=0;p<ploidy;p++){
			if(data->mGenotype[i][snp1][p]!=0)
				data->mGenotype[i][snp2][p]=3-data->mGenotype[i][snp1][p];
		}
	}
}

int main(){
	int sampleNum = 2000;
	int snpNum = 5;
	int ploidy = 2;
	boost::shared_ptr<SHEsisData> data=GenerateRandomData(sampleNum,snpNum,ploidy);
	data->statCount();

//	  std::cout << "Genotype Matrix:\n";
//	  for (int iSample = 0; iSample < sampleNum; iSample++) {
//		  std::cout<<iSample<<"\t"<<(data->vQuantitativeTrait[iSample])<<" ";
//	    for (int iSnp = 0; iSnp < snpNum; iSnp++) {
//	      for (int iChrset = 0; iChrset < ploidy; iChrset++) {
//	        std::cout << data->mGenotype[iSample][iSnp][iChrset] << "/";
//	      }
//	      std::cout << " ";
//	    }
//	    std::cout << "\n";
//	  }

	GeneInteractionQTL gib(data);
	gib.setBinNum(100);
	gib.setSamplePerBin(3);
	gib.setMinBin(2);
	gib.setlb(2);
	gib.setub(3);
	gib.CalGeneInteraction();
	gib.print();
	return 0;
}

