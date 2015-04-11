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
#include <fstream>
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
			if(data->mGenotype[i][snp1][p]!=0)
				data->mGenotype[i][snp2][p]=3-data->mGenotype[i][snp1][p];
		}
//		std::cout<<"sample "<<i<<",after:"<<data->mGenotype[i][snp2][0]<<data->mGenotype[i][snp2][1];

	}
}

bool equal( boost::multi_array<short, 3> data,int sample,int snp1,int snp2, int ploidy){
//	std::cout<<"snp1:";
//	for(int i=0;i<ploidy;i++){
//		std::cout<<data[sample][snp1][i];
//	};
//	std::cout<<"\nsnp2:";
//	for(int i=0;i<ploidy;i++){
//		std::cout<<data[sample][snp2][i];
//	};
//	std::cout<<"\n";
	for(int p=0;p<ploidy;p++){
		if(data[sample][snp1][p]!=data[sample][snp2][p]){
//			std::cout<<"not equal\n";
			return false;
		}
	}
//	std::cout<<"equal\n";
	return true;
}

void SetGXGXG(boost::shared_ptr<SHEsis::SHEsisData> data,int snp1,int snp2,int snp3,double pre){
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
		if(!equal(data->mGenotype,i,snp1,snp2,ploidy)){
//			std::cout<<"sample "<<i<<" set to 1/1\n";
			for(int p=0;p<ploidy;p++){
				data->mGenotype[i][snp3][p]=1;
			}
		}
		else{
//			std::cout<<"sample "<<i<<" set to 2/2\n";
			for(int p=0;p<ploidy;p++){
				data->mGenotype[i][snp3][p]=2;
			}
		}

	}
}

boost::shared_ptr<SHEsis::SHEsisData> GenerateRandomData(int sampleNum,
                                                         int snpNum,
                                                         int chrSetNum) {
  double probabilities[] = {0, 0.5,
                            0.5};  // 0.04 missing phenotype for individuals
  double probabilities2[] = {0.0, 0.5,
                             0.5};  // 0.01 is missing genotype for individuals
  boost::random::discrete_distribution<> dist(probabilities);
  boost::random::discrete_distribution<> dist2(probabilities2);
  boost::shared_ptr<SHEsis::SHEsisData> data(
      new SHEsis::SHEsisData(sampleNum, snpNum, chrSetNum));
  for (int iSample = 0; iSample < sampleNum; iSample++) {
//    BOOST_ASSERT(iSample < data->vLabel.size());
    data->vLabel.push_back(((SHEsis::SampleStatus)dist(boost_rng)));
  }
    for (int iSnp = 0; iSnp < snpNum; iSnp++) {
  	  double maf=dist(boost_rng)==1?0.2:0.4;
  	  double GenoProb[] = {0.0, maf,1-maf};  // 0.01 is missing genotype for individuals
  	  boost::random::discrete_distribution<> distGeno(GenoProb);
  	for (int iSample = 0; iSample < sampleNum; iSample++) {
      for (int iChrset = 0; iChrset < chrSetNum; iChrset++) {
        data->mGenotype[iSample][iSnp][iChrset] = distGeno(boost_rng);
      }
    }
  }
  //	BOOST_CHECK(sampleNum==data->mGenotype.shape()[0]);
  //	BOOST_CHECK(snpNum==data->mGenotype.shape()[1]);
  //	BOOST_CHECK(chrSetNum==data->mGenotype.shape()[2]);
  return data;
}

int main(int argc,char *argv[]){
	string snp[6][3][2]={
			{{"1","2"},{"1","2"},{"1","2"}},
			{{"1","2"},{"1","2"},{"1","1"}},
			{{"1","1"},{"1","1"},{"2","1"}},
			{{"1","2"},{"2","1"},{"1" ,"1"}},
			{{"1","1"},{"2","1"},{"1","1"}},
			{{"1","1"},{"2","2"},{"1","2"}}
	};

	double res[6]={1,1,1,2,2,2};
	int sampleNum = 1000;
	int snpNum = 2000;
	int ploidy = 2;
	sampleNum=atoi(argv[1]);
	snpNum=atoi(argv[2]);
	ploidy=atoi(argv[3]);
	double pre=0.3;
	pre=atof(argv[4]);
	int permutation=atoi(argv[5]);
	double maf=atof(argv[6]);
	int lowb=atoi(argv[7]);
	int hib=atoi(argv[8]);
//	boost::shared_ptr<SHEsisData> data(new SHEsis::SHEsisData(sampleNum, snpNum, ploidy));
	boost::shared_ptr<SHEsisData> data=GenerateRandomData(sampleNum,snpNum,ploidy);
//	for(int i=0;i<snpNum;i=i+2){
//			SetGXG(data,i,i+1,pre);
//	}
	for(int i=0;i<snpNum;i=i+3){
			SetGXGXG(data,i,i+1,i+2,pre);
	}
	data->statCount();
	std::ofstream ped,map;
	std::stringstream name;
	name<<"binary_"<<(snpNum/2)<<"interaction_"<<sampleNum<<"samples_"<<ploidy<<"ploidy_"<<pre<<"pre_";
//	name<<"binary_"<<snpNum<<"snps_"<<sampleNum<<"samples_"<<ploidy<<"ploidy_"<<maf<<"maf";
	std::string pedname=name.str()+".ped";
	std::string mapname=name.str()+".map";
	ped.open(pedname.data());
	map.open(mapname.data());
//	  std::cout << "Genotype Matrix:\n";
	  for (int iSample = 0; iSample < sampleNum; iSample++) {
		  ped<<iSample<<" 0 0 0 1 "<<(data->vLabel[iSample])<<" ";
	    for (int iSnp = 0; iSnp < snpNum; iSnp++) {
	      for (int iChrset = 0; iChrset < ploidy; iChrset++) {
	    	  ped << data->mGenotype[iSample][iSnp][iChrset] << " ";
	      }
	      ped << " ";
	    }
	    ped<< "\n";
	  }
	  for(int snp=0;snp<snpNum;snp++){
		  map<<"1 snp"<<snp<<" 0 "<<snp<<"\n";
	  }
	  ped.close();
	  map.close();
//	SetGXG(data,0,1,0.5);
//	SetGXGXG(data,0,1,2,0.1);
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
	gib.setlb(lowb);
	gib.setPermutation(permutation);
	gib.setub(hib);
	gib.CalGeneInteraction();
	gib.print();
	return 0;
}
