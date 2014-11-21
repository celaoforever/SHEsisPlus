#include "SHEsisData.h"
#include "utility.h"
#include "QTL.h"
#include <iostream>
#include <boost/test/minimal.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/assert.hpp>
#include <boost/random/mersenne_twister.hpp>
boost::mt19937 boost_rng;
boost::shared_ptr<SHEsis::SHEsisData> GenerateRandomData(int sampleNum,
                                                         int snpNum,
                                                         int chrSetNum) {
  double probabilities[] = {0.04, 0.48,
                            0.48};  // 0.04 missing phenotype for individuals
  double probabilities2[] = {0.01, 0.33, 0.33,
                             0.33};  // 0.01 is missing genotype for individuals
                                     //  double probabilities2[] = {0.01, 0.5,
  //                             0.5};  // 0.01 is missing genotype for
  // individuals
  boost::random::uniform_int_distribution<> distQTL(1, 100);
  boost::random::discrete_distribution<> dist(probabilities);
  boost::random::discrete_distribution<> dist2(probabilities2);
  boost::shared_ptr<SHEsis::SHEsisData> data(
      new SHEsis::SHEsisData(sampleNum, snpNum, chrSetNum));
  for (int iSample = 0; iSample < sampleNum; iSample++) {
    BOOST_ASSERT(iSample < data->vLabel.size());
    data->vLabel[iSample] = ((SHEsis::SampleStatus)dist(boost_rng));
    data->vQuantitativeTrait.push_back((double)distQTL(boost_rng) / 20);
    for (int iSnp = 0; iSnp < snpNum; iSnp++) {
      for (int iChrset = 0; iChrset < chrSetNum; iChrset++) {
        data->mGenotype[iSample][iSnp][iChrset] = dist2(boost_rng);
      }
    }
  }
  //	BOOST_CHECK(sampleNum==data.mGenotype.shape()[0]);
  //	BOOST_CHECK(snpNum==data.mGenotype.shape()[1]);
  //	BOOST_CHECK(chrSetNum==data.mGenotype.shape()[2]);
  return data;
}

int test_main(int, char * []) {
  int sampleNum = 400;
  int snpNum = 20;
  int chrSetNum = 2;
  boost::shared_ptr<SHEsis::SHEsisData> testdata =
      GenerateRandomData(sampleNum, snpNum, chrSetNum);
  //  boost::shared_ptr<SHEsis::SHEsisData> testdata(new
  // SHEsis::SHEsisData(6,2,2));
  //  short data[6][2][2]={
  //		  {{1,1},{2,3}},
  //		  {{1,4},{2,3}},
  //		  {{4,4},{3,3}},
  //		  {{1,4},{2,2}},
  //		  {{4,4},{2,3}},
  //		  {{4,4},{2,2}}};
  // for(int sample=0;sample<6;sample++){
  //	 for(int snp=0;snp<2;snp++){
  //		 for(int p=0;p<2;p++){
  //			 testdata->mGenotype[sample][snp][p]=data[sample][snp][p];
  //		 }
  //	 }
  // }
  // double qtl[6]={4,1,1,2,2,2};
  // for(int sample=0;sample<6;sample++){
  //	 testdata->vQuantitativeTrait.push_back(qtl[sample]);
  // }
  std::cout << "Genotype Matrix:\n";
  for (int iSample = 0; iSample < testdata->getSampleNum(); iSample++) {
    for (int iSnp = 0; iSnp < testdata->getSnpNum(); iSnp++) {
      for (int iChrset = 0; iChrset < testdata->getNumOfChrSet(); iChrset++) {
        std::cout << testdata->mGenotype[iSample][iSnp][iChrset] << "/";
      }
      std::cout << " ";
    }
    std::cout << "\n";
  }
  std::cout << "vQuantitativeTrait:\n";
  for (int i = 0; i < testdata->vQuantitativeTrait.size(); i++)
    std::cout << testdata->vQuantitativeTrait[i] << " ";
  std::cout << std::endl;
  //  testdata->statCount();
  //  testdata->printLocusInfo();

  SHEsis::QTL QTLHandle(testdata);
  QTLHandle.setPermutation(1000);
  QTLHandle.QTLPermutation();
  QTLHandle.printRes();
  return boost::exit_success;
}
