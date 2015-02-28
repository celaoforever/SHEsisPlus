/*
 * GeneInteractionQTL_test.cpp
 *
 *  Created on: Jan 10, 2015
 *      Author: ionadmin
 */

#include "GeneInteractionQTL.h"
#include <sstream>
using namespace std;
using namespace SHEsis;
#include <boost/random/discrete_distribution.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <fstream>

boost::mt19937 boost_rng;

bool sortSampleByQtl(const qtl2sampleIdx& v1, const qtl2sampleIdx& v2){
	return (v1.qtl<v2.qtl);
}

boost::shared_ptr<SHEsis::SHEsisData> GenerateRandomData(int sampleNum,
                                                         int snpNum,
                                                         int chrSetNum,double maf) {
  double GenoProb[] = {0.0, maf,1-maf};  // 0.01 is missing genotype for individuals
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
	int binidx=1;
	double unit=1/(double)NumBin/(double)NumBin; //square
//	double unit=1/(double)NumBin; //linear
//	double unit=1/sqrt((double)NumBin);//sqrt
	for(std::list<bin>::iterator iter=bins.begin();iter!=bins.end();iter++){
		double final_unit=unit;
		int interactionNum=0;
		double noise=(double)dice()/100.;
		if(noise<0.0){
			double per=(double)dice()/100.;
			per=per-0.5;
			final_unit=unit*(1+per);
		}
	//	std::cout<<"calc bin "<<binidx<<"\tpre\t"<<((double)binidx*binidx*final_unit)<<"\n";
		for(int i=0;i<iter->SampleIdx.size();i++){
			double n=(double)dice()/100.;
//			double noise=(double)dice()/100;
//			if(noise<0.05){
//				continue;
//			}
			if(n>pre){
				continue;
			}
			n=(double)dice()/100.;
//			if(n<(double)binidx*final_unit) //linear
			if(n<(double)binidx*binidx*final_unit) //square
//			if(n<sqrt((double)binidx)*final_unit) //sqrt
			{
 				continue;
			};

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

int main(int argc,char *argv[]){
	int sampleNum = 1000;
	int snpNum = 2000;
	int ploidy = 3;
	double pre=0.05;
	sampleNum=atoi(argv[1]);
	snpNum=atoi(argv[2]);
	ploidy=atoi(argv[3]);
	pre=atof(argv[4]);
	int permutation=atoi(argv[5]);
	double maf=atof(argv[6]);
	int lowb=2;//atoi(argv[7]);
	int hib=2;//atoi(argv[8]);
	double a=0;
	std::cout<<"sample:"<<sampleNum<<",snp:"<<snpNum<<",ploidy:"<<ploidy<<",pre:"<<pre<<",maf:"<<maf<<",lowb"<<lowb<<",ub"<<hib<<"\n";
	boost::shared_ptr<SHEsisData> data=GenerateRandomData(sampleNum,snpNum,ploidy,maf);
	for(int i=0;i<snpNum;i=i+2){
		std::vector<int> snp;
		snp.push_back(i);
		snp.push_back(i+1);
		SetGXG(data,snp,40,2,pre);
	}
//			std::vector<int> snp;
//			snp.push_back(0);
//			snp.push_back(1);
//			SetGXG(data,snp,40,2,0.5);
	data->statCount();
//	std::ofstream output;
//	output.open("1000samples10snpinteraciton0-1p2.txt");
//		  for (int iSample = 0; iSample < sampleNum; iSample++) {
//			  output<<iSample<<" "<<(data->vQuantitativeTrait[iSample])<<" ";
//		    for (int iSnp = 0; iSnp < snpNum; iSnp++) {
//		      for (int iChrset = 0; iChrset < ploidy; iChrset++) {
//		    	  output << data->mGenotype[iSample][iSnp][iChrset] << " ";
//		      }
//		      output << " ";
//		    }
//		    output<< "\n";
//		  }
	std::ofstream ped,map;
	std::stringstream name;
	name<<"square_"<<(snpNum/2)<<"interaction_"<<sampleNum<<"samples_"<<ploidy<<"ploidy_"<<pre<<"pre_"<<maf<<"maf";
	//name<<"normal_"<<snpNum<<"snps_"<<sampleNum<<"samples_"<<ploidy<<"ploidy_"<<maf<<"maf";
	std::string pedname=name.str()+".ped";
	std::string mapname=name.str()+".map";
	ped.open(pedname.data());
	map.open(mapname.data());
//	  std::cout << "Genotype Matrix:\n";
	  for (int iSample = 0; iSample < sampleNum; iSample++) {
		  ped<<iSample<<" 0 0 0 1 "<<(data->vQuantitativeTrait[iSample])<<" ";
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

	GeneInteractionQTL gib(data);
	gib.setBinNum(2);
//	gib.setSamplePerBin(40);
	gib.setMinBin(0);
	gib.setlb(lowb);
	gib.setPermutation(permutation);
	gib.setub(hib);
	gib.CalGeneInteraction();
	gib.print();
	return 0;
}

