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
#include <boost/assert.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
boost::mt19937 boost_rng;


boost::shared_ptr<SHEsis::SHEsisData> GenerateRandomData(int sampleNum,
                                                         int snpNum,
                                                         int chrSetNum) {
  double probabilities[] = {0, 0.5,
                            0.5};  // 0.04 missing phenotype for individuals
  double GenoProb[] = {0.0, 0.5,
                             0.5};  // 0.01 is missing genotype for individuals

  boost::normal_distribution<> nd(10,2);
  boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > var_nor(rng, nd);
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

int main(){
	int sampleNum = 2000;
	int snpNum = 5;
	int ploidy = 2;
	boost::shared_ptr<SHEsisData> data=GenerateRandomData(sampleNum,snpNum,ploidy);
	data->statCount();

//	  std::cout << "Genotype Matrix:\n";
//	  for (int iSample = 0; iSample < sampleNum; iSample++) {
//		  std::cout<<((int)data->vLabel[iSample]==2?"case: ":"ctrl: ");
//	    for (int iSnp = 0; iSnp < snpNum; iSnp++) {
//	      for (int iChrset = 0; iChrset < ploidy; iChrset++) {
//	        std::cout << data->mGenotype[iSample][iSnp][iChrset] << "/";
//	      }
//	      std::cout << " ";
//	    }
//	    std::cout << "\n";
//	  }

	GeneInteractionQTL gib(data);
	gib.setlb(2);
	gib.setub(3);
	gib.CalGeneInteraction();
	gib.print();
	return 0;
}

