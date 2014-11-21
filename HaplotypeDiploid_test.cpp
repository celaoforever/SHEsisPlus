/*
 * Haplotype_test.cpp
 *
 *  Created on: Aug 19, 2014
 *      Author: ada
 */

#include "HaplotypeDiploid.h"
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
  //	boost::shared_ptr<int> sp(new short[SnpNum]);
  std::vector<boost::shared_ptr<short[]> > res;
  double probabilities[] = {0, 0.3, 0.2, 0.1, 0, 1};  // allele freq
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

boost::shared_ptr<SHEsis::SHEsisData> GenerateHaploData(
    int sampleNum, int snpNum, int chrSetNum,
    std::vector<boost::shared_ptr<short[]> > hap, std::vector<int>& hapcount) {
  boost::shared_ptr<SHEsis::SHEsisData> data(
      new SHEsis::SHEsisData(sampleNum, snpNum, chrSetNum));
  boost::random::uniform_int_distribution<> r(0, hap.size() - 1);
  boost::random::uniform_int_distribution<> r01(0, 1);
  //	std::cout<<"hap.size()="<<hap.size()<<"\n";

  double probabilities[] = {0, 0.48,
                            0.48};  // 0.04 missing phenotype for individuals
  boost::random::discrete_distribution<> dist(probabilities);
  hapcount.resize(hap.size(), 0);
  for (int i = 0; i < sampleNum; i++) {
    data->vLabel[i] = ((SHEsis::SampleStatus)dist(boost_rng));
    int idx1 = r(boost_rng);
    int idx2 = r(boost_rng);
    hapcount[idx1]++;
    hapcount[idx2]++;
    for (int n = 0; n < snpNum; n++) {
      int p = r01(boost_rng);
      data->mGenotype[i][n][p] = hap[idx1][n];
      data->mGenotype[i][n][1 - p] = hap[idx2][n];
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
  int sampleNum = 1000;
  int snpNum = 100;
  int chrSetNum = 2;
  int HapNum = 4;
  std::vector<boost::shared_ptr<short[]> > haps =
      createHaplotype(HapNum, snpNum);
  std::vector<int> hapcount;
  boost::shared_ptr<SHEsis::SHEsisData> data =
      GenerateHaploData(sampleNum, snpNum, chrSetNum, haps, hapcount);
  data->statCount(data->vLabel);
  std::vector<short> mask(3);
  mask[0] = 0;
  mask[1] = 1;
  mask[2] = 1;
  //	data->printLocusInfo();
  //	SHEsis::Haplotype hp(data,2,mask);
  SHEsis::HaplotypeDiploid hp(data);  //,2,mask);
  hp.setSilent(false);
  hp.startHaplotypeAnalysis();
  hp.AssociationTest();
}

int
    // main(){
    test_main(int, char * []) {

  testHp();
  //	return 0;
  return boost::exit_success;
}
