/*
 * SHEsisData_test.cpp
 *
 *  Created on: Aug 8, 2014
 *      Author: ada
 */

#include "SHEsisData.h"
#include "utility.h"
#include <iostream>
#include <boost/test/minimal.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/assert.hpp>
#include <boost/random/mersenne_twister.hpp>
boost::mt19937 boost_rng;
SHEsis::SHEsisData GenerateRandomData(int sampleNum, int snpNum,
                                      int chrSetNum) {
  double probabilities[] = {0.04, 0.48,
                            0.48};  // 0.04 missing phenotype for individuals
  double probabilities2[] = {0.01, 0.33, 0.33,
                             0.33};  // 0.01 is missing genotype for individuals
  boost::random::discrete_distribution<> dist(probabilities);
  boost::random::discrete_distribution<> dist2(probabilities2);
  SHEsis::SHEsisData data(sampleNum, snpNum, chrSetNum);
  for (int iSample = 0; iSample < sampleNum; iSample++) {
    BOOST_ASSERT(iSample < data.vLabel.size());
    data.vLabel[iSample] = ((SHEsis::SampleStatus)dist(boost_rng));
    for (int iSnp = 0; iSnp < snpNum; iSnp++) {
      for (int iChrset = 0; iChrset < chrSetNum; iChrset++) {
        data.mGenotype[iSample][iSnp][iChrset] = dist2(boost_rng);
      }
    }
  }
  //	BOOST_CHECK(sampleNum==data.mGenotype.shape()[0]);
  //	BOOST_CHECK(snpNum==data.mGenotype.shape()[1]);
  //	BOOST_CHECK(chrSetNum==data.mGenotype.shape()[2]);
  return data;
}

int test_main(int, char * []) {
  int sampleNum = 40;
  int snpNum = 2;
  int chrSetNum = 3;
  SHEsis::SHEsisData testdata =
      GenerateRandomData(sampleNum, snpNum, chrSetNum);
  std::cout << "Genotype Matrix:\n";
  for (int iSample = 0; iSample < sampleNum; iSample++) {
    for (int iSnp = 0; iSnp < snpNum; iSnp++) {
      for (int iChrset = 0; iChrset < chrSetNum; iChrset++) {
        std::cout << testdata.mGenotype[iSample][iSnp][iChrset] << "/";
      }
      std::cout << " ";
    }
    std::cout << "\n";
  }
  std::cout << "vLabel:\n";
  for (int i = 0; i < testdata.vLabel.size(); i++)
    std::cout << testdata.vLabel[i] << " ";
  std::cout << std::endl;
  testdata.statCount(testdata.vLabel);
  testdata.printLocusInfo();
  return boost::exit_success;
}
