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
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/shared_ptr.hpp>
boost::mt19937 boost_rng;
std::vector<boost::shared_ptr<short[]> > createHaplotype(int HapNum,
                                                         int SnpNum) {
  std::vector<boost::shared_ptr<short[]> > res;
  double probabilities[] = {0, 0.3, 0.4};  // allele freq
  boost::random::discrete_distribution<> dist(probabilities);
  for (int i = 0; i < HapNum; i++) {
    boost::shared_ptr<short[]> sp(new short[SnpNum]);
    for (int j = 0; j < SnpNum; j++) {
      sp[j] = (short)dist(boost_rng);
    }  //
    res.push_back(sp);
  }

  return res;
}
boost::shared_ptr<SHEsis::SHEsisData> GenerateRandomData(int sampleNum,
                                                         int snpNum,
                                                         int chrSetNum) {
  double probabilities[] = {0, 0.48,
                            0.48};  // 0.04 missing phenotype for individuals
  double probabilities2[] = {0, 0.5,
                             0.5};  // 0.01 is missing genotype for individuals
  boost::random::discrete_distribution<> dist(probabilities);
  boost::random::discrete_distribution<> dist2(probabilities2);
  boost::shared_ptr<SHEsis::SHEsisData> data(
      new SHEsis::SHEsisData(sampleNum, snpNum, chrSetNum));
  for (int iSample = 0; iSample < sampleNum; iSample++) {
    BOOST_ASSERT(iSample < data->vLabel.size());
    data->vLabel[iSample] = ((SHEsis::SampleStatus)dist(boost_rng));
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

boost::shared_ptr<SHEsis::SHEsisData> GenerateHaploData(
    int sampleNum, int snpNum, int chrSetNum,
    std::vector<boost::shared_ptr<short[]> > hap, std::vector<int>& hapcount) {
  boost::shared_ptr<SHEsis::SHEsisData> data(
      new SHEsis::SHEsisData(sampleNum, snpNum, chrSetNum));
  boost::random::uniform_int_distribution<> r(0, hap.size() - 1);

  double probabilities[] = {0, 0.48,
                            0.48};  // 0.04 missing phenotype for individuals
  boost::random::discrete_distribution<> dist(probabilities);
  hapcount.resize(hap.size(), 0);
  for (int i = 0; i < sampleNum; i++) {
    data->vLabel[i] = ((SHEsis::SampleStatus)dist(boost_rng));
    for (int k = 0; k < chrSetNum; k++) {
      int idx = r(boost_rng);
      hapcount[idx]++;
      for (int n = 0; n < snpNum; n++) {
        data->mGenotype[i][n][k] = hap[idx][n];
      }
    }
  }
  std::cout << "haplotype:\n";
  for (int i = 0; i < hap.size(); i++) {
    for (int j = 0; j < snpNum; j++) {
      std::cout << hap[i][j];
    };
    std::cout << "\t" << hapcount[i] << "\n";
  };
  return data;
}

void testHp() {
  int sampleNum = 1200;
  int snpNum = 2;
  int chrSetNum = 2;
  int HapNum = 3;
  std::vector<boost::shared_ptr<short[]> > haps =
      createHaplotype(HapNum, snpNum);
  std::vector<int> hapcount;
  boost::shared_ptr<SHEsis::SHEsisData> data =
      GenerateHaploData(sampleNum, snpNum, chrSetNum, haps, hapcount);
  //  data->statCount(data->vLabel);
  SHEsis::Haplotype hp(data);
  hp.setSilent(false);
  hp.startHaplotypeAnalysis();
  hp.AssociationTest();
}

void testAs() {
  boost::shared_ptr<int[]> sizes(new int[3]);
  sizes[0] = 2;
  sizes[1] = 1;
  sizes[2] = 2;

  boost::shared_ptr<int[]> array(new int[2]);
  array[0] = 1;
  array[1] = 1;

  SHEsis::ArrayStorage AS1(sizes);
  SHEsis::ArrayStorage AS2(sizes, array);
  AS2.getArray();
  AS2.getDimension();
  AS2.getSizes();
  AS2.set(1, 1);
}

void testIV() {
  SHEsis::IndexingVariables iv;
  iv.add("test", SetSharedPtr(2, 1, 1));
  iv.getEnumeration("test", SetSharedPtr(1, 0));
  iv.printHmKey();
}
int test_main(int, char * []) {

  testHp();
  return boost::exit_success;
}
