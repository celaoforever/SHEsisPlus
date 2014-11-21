/*
 * HWETest_test.cpp
 *
 *  Created on: Aug 11, 2014
 *      Author: ada
 */
#include "HWETest.h"
#include "utility.h"
#include <iostream>
#include <boost/test/minimal.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/assert.hpp>
#include <boost/random/mersenne_twister.hpp>
boost::mt19937 boost_rng;

boost::shared_ptr<SHEsis::SHEsisData> GenerateRandomData(int sampleNum,
                                                         int snpNum,
                                                         int chrSetNum) {
  double probabilities[] = {0, 0.5,
                            0.5};  // 0.04 missing phenotype for individuals
  double probabilities2[] = {0, 0.5,
                             0.5};  // 0.01 is missing genotype for individuals
  boost::random::discrete_distribution<> dist(probabilities);
  boost::random::discrete_distribution<> dist2(probabilities2);
  boost::shared_ptr<SHEsis::SHEsisData> data(
      new SHEsis::SHEsisData(sampleNum, snpNum, chrSetNum));
  for (int iSample = 0; iSample < sampleNum; iSample++) {
    BOOST_ASSERT(iSample < data->vLabel.size());
    data->vLabel[iSample] = ((SHEsis::SampleStatus)dist(boost_rng));
    //    data->vQuantitativeTrait.push_back(dist(boost_rng));
    for (int iSnp = 0; iSnp < snpNum; iSnp++) {
      for (int iChrset = 0; iChrset < chrSetNum; iChrset++) {
        data->mGenotype[iSample][iSnp][iChrset] = dist2(boost_rng);
      }
    }
  }
  //	BOOST_CHECK(sampleNum==data->mGenotype.shape()[0]);
  //	BOOST_CHECK(snpNum==data->mGenotype.shape()[1]);
  //	BOOST_CHECK(chrSetNum==data->mGenotype.shape()[2]);
  return data;
}

int test_main(int, char * []) {
  int sampleNum = 20000;
  int snpNum = 1;
  int chrSetNum = 2;
  boost::shared_ptr<SHEsis::SHEsisData> testdata =
      GenerateRandomData(sampleNum, snpNum, chrSetNum);
  std::cout << "Genotype Matrix:\n";
  for (int iSample = 0; iSample < sampleNum; iSample++) {
    for (int iSnp = 0; iSnp < snpNum; iSnp++) {
      for (int iChrset = 0; iChrset < chrSetNum; iChrset++) {
        std::cout << testdata->mGenotype[iSample][iSnp][iChrset] << "/";
      }
      std::cout << " ";
    }
    std::cout << "\n";
  }
  //  std::cout << "vLabel:\n";
  //  for (int i = 0; i < testdata->vLabel.size(); i++)
  //    std::cout << testdata->vLabel[i] << " ";
  ////  testdata->vQuantitativeTrait.push_back(1);
  //  std::cout << std::endl;

  SHEsis::HWETest HWETestHandler(testdata);
  // HWETestHandler.data->statCount(HWETestHandler.data->vLabel);
  //  HWETestHandler.data->statCount();
  //  HWETestHandler.data->printLocusInfo();
  HWETestHandler.AllSnpHWETest();
  HWETestHandler.printHWETestResults();

  return boost::exit_success;
}
