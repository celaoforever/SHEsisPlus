/*
 * GeneInteraction.cpp
 *
 *  Created on: Dec 30, 2014
 *      Author: ionadmin
 */

#include "GeneInteraction.h"
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <boost/assert.hpp>
namespace SHEsis {


bool GeneInteraction::genotypeEqual(std::string geno,int sample,int snp){
	 return (this->mGenotypeStr[sample][snp]==geno);
}

GeneInteraction::GeneInteraction(boost::shared_ptr<SHEsisData> data):lowbound(2),upperbound(2),data(data),
		mGenotypeStr(boost::extents[data->getSampleNum()][data->getSnpNum()]){
	if(this->data->vLocusInfo[0].BothGenotypeCount.begin() == this->data->vLocusInfo[0].BothGenotypeCount.end()){
		this->data->statCount(this->data->vLabel);
	}
	for(int sample=0;sample<this->data->getSampleNum();sample++){
		for(int snp=0;snp<this->data->getSnpNum();snp++){
			std::stringstream geno("");
			geno<<(this->data->mGenotype[sample][snp][0]);
			for(int p=1;p<this->data->getNumOfChrSet();p++){
				geno<<"/"<<this->data->mGenotype[sample][snp][p];
			}
			this->mGenotypeStr[sample][snp]=geno.str();
		}
	}
}

void GeneInteraction::GenerateSNPCombination(int snpnum,std::vector<int> snpidx,std::vector<std::vector<int> >& ret){
	ret.clear();
	std::vector<int> r;
//    std::cout<<"GenerateSNPCombination()\nsnpnum="<<snpnum<<"\n";
//    std::cout<<"snpidx:";
//    for(int i=0;i<snpidx.size();i++)
//    	std::cout<<snpidx[i]<<",";
//    std::cout<<"\n";
	BOOST_ASSERT(snpnum<=snpidx.size());
	for(int i=0; i <snpnum;i++)
		r.push_back(snpidx[i]);
	  do {
	    std::vector<int> tmp = r;
	    ret.push_back(tmp);
	  } while (next_combination(snpidx.begin(), snpidx.end(), r.begin(), r.end()));
}

void GeneInteraction::GenerateSNPCombination(int snpnum,std::vector<std::vector<int> >& ret){
	ret.clear();
	std::vector<int> r;
	std::vector<int> c;

	BOOST_ASSERT(snpnum<=this->data->getSnpNum());
//	BOOST_ASSERT(snpnum<=total);
	for(int i=0; i <snpnum;i++)
		r.push_back(i);
//	std::cout<<"snpnum:"<<this->data->getSnpNum()<<"\n";
	for(int i=0;i<this->data->getSnpNum();i++)
//	for(int i=0;i<total;i++)
		c.push_back(i);

	  do {
	    std::vector<int> tmp = r;
	    ret.push_back(tmp);
	  } while (next_combination(c.begin(), c.end(), r.begin(), r.end()));
}

double GeneInteraction::getInformationInteraction(std::vector<int> samples,std::vector<int> snps){
	double ret=0;
//	std::cout<<"getInformationInteraction()\nsnps:";
//	for(int i=0;i<snps.size();i++){
//		std::cout<<snps[i]<<",";
//	}
	for(int i=1;i<=snps.size();i++){
//		std::cout<<"number of snp:"<<i<<"\n";
		std::vector<std::vector<int> > combination;
		this->GenerateSNPCombination(i,snps,combination);
		for(int j =0;j<combination.size();j++){
			ret+=pow(-1.0,(double)(snps.size()-i))*(this->getEntropy(samples,combination[j]));
		}
	}
	return ret;
}

double GeneInteraction::getEntropy(std::vector<int> samples, std::vector<int> snps){
	double ret=0;
//	std::cout<<"getEntropy()\nsnps:";
//	for(int i=0;i<snps.size();i++){
//		std::cout<<snps[i]<<",";
//	}
	std::vector<std::vector<std::string> > cp;
	this->GenerateGenotypeCombination(snps,cp);

//	std::cout<<"\ncp:\n";
//	for(int i=0;i<cp.size();i++){
//		for(int j=0;j<cp[i].size();j++){
//			std::cout<<cp[i][j]<<",";
//		}
//		std::cout<<"\n";
//	}
	for(int cpIdx=0;cpIdx<cp.size();cpIdx++){
		int count=0;
		for(int sample=0;sample<samples.size();sample++){
			bool equal=true;
			for(int i=0;i<cp[cpIdx].size();i++){
				if(!this->genotypeEqual(cp[cpIdx][i],samples[sample],snps[i])){
					equal = false;
					break;
				}
			}
			if(equal)
				count++;
		}
//		std::cout<<"genotype combination for case:\n";
//		for(int i=0;i<cp[cpIdx].size();i++){
//			std::cout<<cp[cpIdx][i]<<",";
//		}

		double rate=(double)count/(double)samples.size();
//		std::cout<<"\nrate="<<rate<<"\n";
		if(rate!=0)
			ret+=rate*log(rate);
	}
	return ret*(-1);
}

void GeneInteraction::GenerateGenotypeCombination(std::vector<int> &Snp,std::vector<std::vector<std::string> >& ret){
	std::vector<std::vector<std::string> > idx;
	ret.clear();
	for(int i=0;i<Snp.size();i++){
		std::vector<std::string> _idx;
		boost::unordered_map<std::string, double>::iterator iter;
		for(iter=this->data->vLocusInfo[Snp[i]].BothGenotypeCount.begin();iter!=this->data->vLocusInfo[Snp[i]].BothGenotypeCount.end();iter++){
			_idx.push_back(iter->first);
		}
		idx.push_back(_idx);
	}
	cart_product(ret,idx);
}

GeneInteraction::~GeneInteraction() {
	// TODO Auto-generated destructor stub
}

} /* namespace SHEsis */
