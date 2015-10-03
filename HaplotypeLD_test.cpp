/*
 * HaplotypeLD_test.cpp
 *
 *  Created on: May 22, 2015
 *      Author: ionadmin
 */

#include "HaplotypeLD.h"
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
  double probabilities[] = {0, 0.7, 0.3};  // allele freq
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

  double probabilities[] = {0.02, 0.49,
                            0.49};  // 0.04 missing phenotype for individuals
  boost::random::discrete_distribution<> dist(probabilities);
  hapcount.resize(hap.size(), 0);
  for (int i = 0; i < sampleNum; i++) {
    data->vLabel.push_back((SHEsis::SampleStatus)dist(boost_rng));
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
  int sampleNum = 1000;
  int snpNum = 10;
  int chrSetNum = 2;
  int HapNum = 10;
  std::vector<boost::shared_ptr<short[]> > haps =
      createHaplotype(HapNum, snpNum);
  std::vector<int> hapcount;
  boost::shared_ptr<SHEsis::SHEsisData> data =
      GenerateHaploData(sampleNum, snpNum, chrSetNum, haps, hapcount);
  data->statCount(data->vLabel);
   /* std::cout<<"Genotype matrix:\n";
    for(int i=0;i<data->getSampleNum();i++){
  	  for(int j=0;j< data->getSnpNum();j++){
  		  for(int p=0;p<data->getNumOfChrSet();p++){
  			  std::cout<<data->mGenotype[i][j][p]<<"/";
  		  }
  		  std::cout<<" ";
  	  }
  	  std::cout<<"\n";
    }
    */
  std::vector<short> mask(5);
  mask[0] = 1;
  mask[1] = 1;
  mask[2] = 1;
  mask[3] = 1;
  mask[4] = 1;
  //	data->printLocusInfo();
  //	SHEsis::Haplotype hp(data,2,mask);
  SHEsis::HaplotypeLD hp(data);  //,2,mask);
  hp.phaseAll();
  hp.AssociationTest();
  hp.setFreqThreshold(0);
  std::cout<<hp.reporttxt();
}

int
    // main(){
    test_main(int, char * []) {

  testHp();
  //	return 0;
  return boost::exit_success;
}
