/*
 * GeneInteractionBinary_test.cpp
 *
 *  Created on: Jan 1, 2015
 *      Author: ionadmin
 */

#include "GeneInteractionBinary.h"

using namespace std;
using namespace SHEsis;
#include <boost/random/discrete_distribution.hpp>
#include <boost/assert.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
boost::mt19937 boost_rng;

void SetGXG(boost::shared_ptr<SHEsis::SHEsisData> data,int snp1,int snp2,double pre){
	int ploidy=data->getNumOfChrSet();
	boost::uniform_int<> dist(1,100);
	boost::variate_generator<boost::mt19937,boost::uniform_int<> > dice(boost_rng,dist);
	for(int i=0;i<data->getSampleNum();i++){
		if(data->vLabel[i]==CONTROL)
			continue;

		double n=(double)dice()/100.;
		if(n>pre){
//			std::cout<<"n="<<n<<",skip sample"<<i<<"\n";
			continue;
		}
//		std::cout<<"sample "<<i<<",origin:"<<data->mGenotype[i][snp2][0]<<data->mGenotype[i][snp2][1];
		for(int p=0;p<ploidy;p++){
			data->mGenotype[i][snp2][p]=3-data->mGenotype[i][snp1][p];
		}
//		std::cout<<"sample "<<i<<",after:"<<data->mGenotype[i][snp2][0]<<data->mGenotype[i][snp2][1];

	}
}

boost::shared_ptr<SHEsis::SHEsisData> GenerateRandomData(int sampleNum,
                                                         int snpNum,
                                                         int chrSetNum) {
  double probabilities[] = {0, 0.5,
                            0.5};  // 0.04 missing phenotype for individuals
  double probabilities2[] = {0.1, 0.45,
                             0.45};  // 0.01 is missing genotype for individuals
  boost::random::discrete_distribution<> dist(probabilities);
  boost::random::discrete_distribution<> dist2(probabilities2);
  boost::shared_ptr<SHEsis::SHEsisData> data(
      new SHEsis::SHEsisData(sampleNum, snpNum, chrSetNum));
  for (int iSample = 0; iSample < sampleNum; iSample++) {
//    BOOST_ASSERT(iSample < data->vLabel.size());
    data->vLabel.push_back(((SHEsis::SampleStatus)dist(boost_rng)));
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

int main(){
	string snp[6][3][2]={
			{{"1","2"},{"1","2"},{"1","2"}},
			{{"1","2"},{"1","2"},{"1","1"}},
			{{"1","1"},{"1","1"},{"2","1"}},
			{{"1","2"},{"2","1"},{"1" ,"1"}},
			{{"1","1"},{"2","1"},{"1","1"}},
			{{"1","1"},{"2","2"},{"1","2"}}
	};

//	double res[6]={1,1,1,2,2,2};
	int sampleNum = 2000;
	int snpNum = 30;
	int ploidy = 2;
//	boost::shared_ptr<SHEsisData> data(new SHEsis::SHEsisData(sampleNum, snpNum, ploidy));
	boost::shared_ptr<SHEsisData> data=GenerateRandomData(sampleNum,snpNum,ploidy);
	SetGXG(data,28,29,0.2);
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
//	genotype and phenotype
//	for(int i=0;i<sampleNum;i++){
//		data->vLabel.push_back((SampleStatus)((int)res[i]));
//		for(int j=0;j<snpNum;j++){
//			for(int p=0;p<ploidy;p++){
//				data->mGenotype[i][j][p]=data->GetAlleleCode(snp[i][j][p]);
//			}
//		}
//	}
//	std::cout<<"snpnum="<<data->getSnpNum()<<"\n";
	GeneInteractionBinary gib(data);
	gib.setlb(3);
	gib.setub(3);
	gib.CalGeneInteraction();
	gib.print();
	return 0;
}
