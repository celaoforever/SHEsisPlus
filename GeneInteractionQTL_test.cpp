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

bool sortSampleByQtl(const qtl2sampleIdx& v1, const qtl2sampleIdx& v2){
	return (v1.qtl<v2.qtl);
}

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

bool equal( boost::multi_array<short, 3> data,int sample,int snp1,int snp2, int ploidy){
	for(int p=0;p<ploidy;p++){
		if(data[sample][snp1][p]!=data[sample][snp2][p]){
			return false;
		}
	}
	return true;
}

void SetGXG(boost::shared_ptr<SHEsis::SHEsisData> data,std::vector<int> snp,int NumBin, int MinSamplesPerBin,double pre){
	int ploidy=data->getNumOfChrSet();
	boost::uniform_int<> dist(1,100);
	boost::variate_generator<boost::mt19937,boost::uniform_int<> > dice(boost_rng,dist);
	//split into bins
	std::list<bin> bins;
	std::vector<qtl2sampleIdx> validSamples;
	for(int i=0;i<data->getSampleNum();i++){
		qtl2sampleIdx s;
		s.idx=i;
		s.qtl=data->vQuantitativeTrait[i];
		validSamples.push_back(s);
	};
	std::sort(validSamples.begin(),validSamples.end(),sortSampleByQtl);
	int maxbin=validSamples.size()/MinSamplesPerBin;
	int SamplesPerBin=validSamples.size()/NumBin;
	if(NumBin>maxbin){
		NumBin=maxbin;
		SamplesPerBin=MinSamplesPerBin;
	}else{
		int rest=validSamples.size()%NumBin;
		if(rest>MinSamplesPerBin)
			NumBin++;
	}

	for(int i=0;i<NumBin-1;i++){
//		std::cout<<"bin "<<i<<":";
		bin b;
		b.meanqtl=0;
		for(int j=0;j<SamplesPerBin;j++){
			int sampleidx=validSamples[i*SamplesPerBin+j].idx;
//			std::cout<<sampleidx<<",";
			b.SampleIdx.push_back(sampleidx);
		}
		bins.push_back(b);
//		std::cout<<"\n";
	}

	//deal with the rest
	bin b;
	for(int j=(NumBin-1)*SamplesPerBin;j<validSamples.size();j++){
		b.SampleIdx.push_back(validSamples[j].idx);
	}
	bins.push_back(b);
	int binidx=0;
	double unit=1/(double)NumBin;
	for(std::list<bin>::iterator iter=bins.begin();iter!=bins.end();iter++){
		int interactionNum=0;
		for(int i=0;i<iter->SampleIdx.size();i++){
			double n=(double)dice()/100.;
			if(n>pre){
				continue;
			}
			n=(double)dice()/100.;
			if(n<(double)binidx*unit)
				continue;

			int sampleidx=iter->SampleIdx[i];
			if(snp.size()== 2 ){

//				std::cout<<"set gene interaction for bin "<<binidx<<",sample "<<sampleidx<<", qtl="<<data->vQuantitativeTrait[sampleidx]<<"\n";
				for(int p=0;p<ploidy;p++){

					if(data->mGenotype[sampleidx][snp[1]][p]!=0)
						interactionNum++;
						data->mGenotype[sampleidx][snp[0]][p]=data->mGenotype[sampleidx][snp[1]][p];
				}
			}else if(snp.size() == 3 ){
				if(!equal(data->mGenotype,sampleidx,snp[0],snp[1],ploidy)){
						for(int p=0;p<ploidy;p++){
						data->mGenotype[sampleidx][snp[2]][p]=1;
					}
				}
				else{
					for(int p=0;p<ploidy;p++){
						data->mGenotype[sampleidx][snp[2]][p]=2;
					}
				}
			}
		}
//		cout<<"interaction num for bin "<<binidx<<":"<<interactionNum<<"\n";
		binidx++;
	}
}

int main(){
	int sampleNum = 1000;
	int snpNum = 4;
	int ploidy = 2;
	boost::shared_ptr<SHEsisData> data=GenerateRandomData(sampleNum,snpNum,ploidy);
	std::vector<int> snp;
	snp.push_back(0);
	snp.push_back(1);
	snp.push_back(2);
	SetGXG(data,snp,3,2,1);
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
	gib.setBinNum(20);
	gib.setSamplePerBin(2);
	gib.setMinBin(0);
	gib.setlb(2);
	gib.setub(3);
	gib.CalGeneInteraction();
	gib.print();
	return 0;
}

