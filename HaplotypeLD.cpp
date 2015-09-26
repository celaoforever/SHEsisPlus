/*
 * HaplotypeLD.cpp
 *
 *  Created on: May 22, 2015
 *      Author: ionadmin
 */

#include "HaplotypeLD.h"
#include "LDTest.h"
#include <boost/algorithm/string.hpp>
namespace SHEsis {

HaplotypeLD::HaplotypeLD(boost::shared_ptr<SHEsisData> data):data(data),curCode(1),ldT(0.1),adjust(false),silent(true),lft(0.03){
	  for (int i = 0; i < this->data->getSnpNum(); i++) {
	    this->SnpIdx.push_back(i);
	  };
};

HaplotypeLD::HaplotypeLD(boost::shared_ptr<SHEsisData> data, int Snp,
              std::vector<short> mask):data(data),curCode(1),ldT(0.1),adjust(false),silent(true),lft(0.03){
	  for (int i = 0; i < mask.size(); i++) {
	    if (mask[i]) {
	      this->SnpIdx.push_back(i);
	    };
	  };
	  this->mask=mask;
};

void HaplotypeLD::getAdjacentLD(){
	  LDTest ld(this->data);
	  if (this->data->vLocusInfo[0].BothAlleleCount.size() == 0 &&
	      this->data->vQuantitativeTrait.size() > 0) {
	    this->data->statCount();
	  } else if (this->data->vLocusInfo[0].BothAlleleCount.size() == 0 &&
	             this->data->vQuantitativeTrait.size() == 0)
	    this->data->statCount(this->data->vLabel);
	  this->ldPattern.reset(new double[this->SnpIdx.size()-1]);
	  for(int i=0;i<this->SnpIdx.size()-1;i++){
		  double R2,D;
		  ld.TwoLociLDTest(SnpIdx[i], SnpIdx[i+1], LD_IN_BOTH, R2, D);
		  this->ldPattern[i]=R2;
	  }
//	  std::cout<<"ld pattern:";
//	  for(int i=0;i<this->SnpIdx.size()-1;i++){
//		  std::cout<<this->ldPattern[i]<<",";
//	  }
//	  std::cout<<"\n";
}

std::vector<int> HaplotypeLD::getLDBlock(double t){
	this->getAdjacentLD();
	std::vector<int> ret;
	ret.push_back(0);
	for(int i=0;i<this->SnpIdx.size()-1;i++){
		if(this->ldPattern[i]<t){
			ret.push_back(i+1);
		}
	}
	ret.push_back(this->SnpIdx.size());
	//snp from SnpIdx[ret[i]] and SnpIdx[ret[i+1]-1] is a block
	this->blockData.reset(new SHEsisData(this->data->getSampleNum(),ret.size()-1,this->data->getNumOfChrSet()));
	return ret;
}

void HaplotypeLD::reducePhasedHap(boost::multi_array<short, 3> PhasedData,int curBlockIdx){
	for(int sample=0;sample<PhasedData.shape()[0];sample++){
		for(int ploidy=0;ploidy<PhasedData.shape()[2];ploidy++){
			std::stringstream ss;
			for(int snp=0;snp<PhasedData.shape()[1];snp++){
				ss<<PhasedData[sample][snp][ploidy]<<",";
			}
			hapMap::map_by<haplo>::const_iterator iter=this->hapcode.by<haplo>().find(ss.str());
			if(this->hapcode.by<haplo>().end()==iter){
				this->hapcode.insert(hapMap::value_type(ss.str(),this->curCode));
				this->blockData->mGenotype[sample][curBlockIdx][ploidy]=this->curCode;
				this->curCode++;
			}else{
				this->blockData->mGenotype[sample][curBlockIdx][ploidy]=iter->get<code>();
			}
		}
	}
}

void HaplotypeLD::phaseAll(){
	std::vector<int> block=this->getLDBlock(this->ldT);
//	std::cout<<"block:";
//	for(int i=0;i<block.size();i++){
//		std::cout<<block[i]<<",";
//	}
//	std::cout<<"\n";
	for(int i=0;i<block.size()-1;i++){
//		std::cout<<"current phased snp"<<this->SnpIdx[block[i]]<<" to "<<this->SnpIdx[block[i+1]-1]<<"\n";
		this->getHaplotypeSub(block[i],block[i+1]-1/*this->SnpIdx[block[i]],this->SnpIdx[block[i+1]-1]*/);
		this->reducePhasedHap(this->hap->PhasedData,i);
	}
	this->blockData->vLabel=this->data->vLabel;
	this->blockData->vQuantitativeTrait=this->data->vQuantitativeTrait;
	this->hap.reset(new HaplotypeEM(this->blockData));
	this->hap->setSilent(this->silent);
	this->hap->setFreqThreshold(this->lft);
	this->hap->setAdjust(this->adjust);
	this->hap->startHaplotypeAnalysis();
	for(int i=0;i<this->hap->Results.haplotypes.size();i++){
		std::string ss;
		boost::shared_ptr<short[]> cur=this->hap->Results.haplotypes[i];
		boost::shared_ptr<short[]> res(new short[this->SnpIdx.size()]);
		int idx=0;
		for(int j=0;j<block.size()-1;j++){
			hapMap::map_by<code>::const_iterator iter=this->hapcode.by<code>().find(cur[j]);
			BOOST_ASSERT(iter!=this->hapcode.by<code>().end());
			ss=iter->get<haplo>();
			std::vector<std::string> vec;
			boost::algorithm::split(vec,ss,boost::algorithm::is_any_of(","));
			for(int k=0;k<vec.size()-1;k++){
				res[idx++]=(short)atoi(vec[k].data());
			}
		}
		this->haplotypes.push_back(res);
	}
    for(int i=0;i<this->data->getSampleNum();i++){
  	  	 std::cerr<<i<<"\t";
  	  for(int j=0;j<this->data->getNumOfChrSet();j++){
  		  		for(int k=0;k<this->SnpIdx.size();k++){
  		  	    	std::cerr<<this->data->getallele(this->haplotypes[this->hap->Results.genotypes[i][j]][k]);
  		  		}
  		  		std::cerr<<"\t";
  		 // std::cout<<this->Results.genotypes[i][j]<<"\t";
  	  }
  	  std::cerr<<"\n";
    }


	//std::cout<<"haplotypes:\n";
	//for(int i=0;i<this->haplotypes.size();i++){
//		for(int j=0;j<this->SnpIdx.size();j++){
//			std::cout<<this->data->getallele(this->haplotypes[i][j]);
//		}
//		std::cout<<":"<< (this->hap->Results.BothCount[i])<<"\n";
//	}

}

void HaplotypeLD::AssociationTest(){
	this->hap->AssociationTest();

}
std::string HaplotypeLD::reporthtml(){
	this->hap->Results.haplotypes=this->haplotypes;
	this->hap->SnpIdx=this->SnpIdx;
	this->hap->data=this->data;
	return this->hap->reporthtml();
}
std::string HaplotypeLD::reporttxt(){
	return this->hap->reporttxt();
}

void HaplotypeLD::getHaplotypeSub(int start, int end){
	std::vector<short> m(this->data->getSnpNum(),0);
	for(int i=start;i<=end;i++){
		m[this->SnpIdx[i]]=1;
	};
	this->hap.reset(new HaplotypeEM(this->data,end-start+1,m));
	this->hap->ShowResults(false);
	this->hap->startHaplotypeAnalysis();
}


HaplotypeLD::~HaplotypeLD() {
	// TODO Auto-generated destructor stub
}

} /* namespace SHEsis */
