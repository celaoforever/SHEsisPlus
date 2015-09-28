/*
 * HaplotypeLD.cpp
 *
 *  Created on: May 22, 2015
 *      Author: ionadmin
 */

#include "HaplotypeLD.h"
#include "LDTest.h"
#include <boost/algorithm/string.hpp>
#define MIN(a, b) a <= b ? a : b
namespace SHEsis {
void printArray(boost::multi_array<short,3> m){
	for(int i=0;i<m.shape()[0];i++){
		for(int j=0;j<m.shape()[1];j++){
			for(int p=0;p<m.shape()[2];p++){
				std::cout<<m[i][j][p];
			}
			std::cout<<" ";
		}
		std::cout<<"\n";
	}
}

void printArray(boost::multi_array<std::string,2> m){
	for(int i=0;i<m.shape()[0];i++){
		for(int j=0;j<m.shape()[1];j++){
			std::cout<<m[i][j]<<"\n";
		}
		std::cout<<"\n";
	}
}
HaplotypeLD::HaplotypeLD(boost::shared_ptr<SHEsisData> data):data(data),numSnpPerBlock(5),curCode(1),ldT(0.1),adjust(false),silent(true),lft(0.03){
	  for (int i = 0; i < this->data->getSnpNum(); i++) {
	    this->SnpIdx.push_back(i);
	  };
};

HaplotypeLD::HaplotypeLD(boost::shared_ptr<SHEsisData> data, int Snp,
              std::vector<short> mask):data(data),numSnpPerBlock(5),curCode(1),ldT(0.1),adjust(false),silent(true),lft(0.03){
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

std::string getPartialHap(short query){

		std::stringstream ss;
		ss<<query<<",";
		return ss.str();
}

std::string getPartialHap(short query, hapMap hapcode){
	hapMap::map_by<code>::const_iterator iter=hapcode.by<code>().find(query);
    if(iter!=hapcode.by<code>().end()){
		return iter->get<haplo>();
	}else
	{
		std::stringstream ss;
		ss<<query<<",";
		return ss.str();
	}
}

void HaplotypeLD::reducePhasedHap(hapMap& hc,boost::shared_ptr<SHEsisData> phased,boost::multi_array<short, 3> PhasedData,int curBlockIdx){
	boost::multi_array<std::string,2> strPhasedData(boost::extents[PhasedData.shape()[0]][PhasedData.shape()[2]]);
	for(int sample=0;sample<PhasedData.shape()[0];sample++){
		for(int ploidy=0;ploidy<PhasedData.shape()[2];ploidy++){
			for(int snp=0;snp<PhasedData.shape()[1];snp++){
				strPhasedData[sample][ploidy]+=getPartialHap(PhasedData[sample][snp][ploidy]);
			}
		}
	}
	//std::cout<<"str Phased Data:\n";
	//printArray(strPhasedData);
	for(int sample=0;sample<PhasedData.shape()[0];sample++){
		for(int ploidy=0;ploidy<PhasedData.shape()[2];ploidy++){
			hapMap::map_by<haplo>::const_iterator iter=hc.by<haplo>().find(strPhasedData[sample][ploidy]);
			if(hc.by<haplo>().end()==iter){
				hc.insert(hapMap::value_type(strPhasedData[sample][ploidy],this->curCode));
				phased->mGenotype[sample][curBlockIdx][ploidy]=this->curCode;
				//std::cout<<strPhasedData[sample][ploidy]<<":"<<this->curCode<<"\n";
				this->curCode++;
			}else{
				phased->mGenotype[sample][curBlockIdx][ploidy]=iter->get<code>();
			}
		}
	}

}

void HaplotypeLD::updateHapCode(hapMap &source, hapMap& output){

	hapMap::map_by<code>::const_iterator iter=output.by<code>().begin();
	hapMap::map_by<code>::const_iterator iter_=source.by<code>().begin();
	/*std::cout<<"source:\n";
	for(;iter_!=source.by<code>().end();iter_++){
		std::cout<<"code:"<<iter_->get<code>()<<",hap:"<<iter_->get<haplo>()<<"\n";
	};
	std::cout<<"dest before:\n";
	for(;iter!=output.by<code>().end();iter++){
		std::cout<<"code:"<<iter->get<code>()<<",hap:"<<iter->get<haplo>()<<"\n";
	};
  */
	iter=output.by<code>().begin();
	for(;iter!=output.by<code>().end();iter++){
		std::vector<std::string> vec;
		std::stringstream ss;
		//std::cout<<"current haplo:"<<iter->get<haplo>()<<"\n";
		boost::algorithm::split(vec,iter->get<haplo>(),boost::algorithm::is_any_of(","));
		for(int i=0;i<vec.size()-1;i++){
			//std::cout<<"current code:"<<(short)atoi(vec[i].c_str())<<"\n";
			hapMap::map_by<code>::const_iterator iter2=source.by<code>().find((short)atoi(vec[i].c_str()));
			BOOST_ASSERT(iter2!=source.by<code>().end());
			ss<<iter2->get<haplo>();
		}
		short c=iter->get<code>();
		std::string hap=iter->get<haplo>();
		output.by<code>().replace_data(output.by<code>().find(c),ss.str());
		//output.by<code>().erase(c);
		//output.insert(hapMap::value_type(ss.str(),c));
	}
	/*std::cout<<"after:\n";
	iter=output.by<code>().begin();
	for(;iter!=output.by<code>().end();iter++){
		std::cout<<"code:"<<iter->get<code>()<<",hap:"<<iter->get<haplo>()<<"\n";
	};
	*/
	source.clear();
}
void HaplotypeLD::phaseAll(){

	this->DataToPhase.reset(new SHEsisData(this->data->getSampleNum(),this->data->getSnpNum(),this->data->getNumOfChrSet()));
	this->DataToPhase->mGenotype.resize(extension(this->data->mGenotype));
	this->DataToPhase->mGenotype=this->data->mGenotype;
	this->DataToPhase->vLabel=this->data->vLabel;
	this->DataToPhase->vQuantitativeTrait=this->data->vQuantitativeTrait;
	int startsnp,endsnp;
	int count=0;
	while(true)
	{
		//std::cout<<"DataToPhase:\n";
		//printArray(this->DataToPhase->mGenotype);
		//std::cout<<"num snp:"<<this->DataToPhase->getSnpNum()<<",per block:"<<this->numSnpPerBlock<<"\n";
		int numOfBlock=(this->DataToPhase->getSnpNum()-1)/this->numSnpPerBlock+1;
		this->PhasedData.reset(new SHEsisData(this->DataToPhase->getSampleNum(),numOfBlock,this->DataToPhase->getNumOfChrSet()));
		this->PhasedData->vLabel=this->data->vLabel;
		this->PhasedData->vQuantitativeTrait=this->data->vQuantitativeTrait;
		for(int i=0;i<numOfBlock;i++){
			startsnp=i*this->numSnpPerBlock;
			endsnp=MIN((1+i)*this->numSnpPerBlock-1,this->DataToPhase->getSnpNum()-1);
			this->getHaplotypeSub(startsnp,endsnp,this->DataToPhase);
			//std::cout<<"Hap PhasedData:\n";
			//printArray(this->hap->PhasedData);
			//std::cout<<"count="<<count<<","<<(count%2?"reduced":"hapcode")<<"\n";
			this->reducePhasedHap(count%2?this->reduced:this->hapcode,this->PhasedData,this->hap->PhasedData,i);
		}
		this->curCode=1;
		if (!this->hapcode.empty()&&!this->reduced.empty()){
			//std::cout<<"count="<<count<<",source="<<(count%2?"hapcode":"reduced")<<"\n";
			this->updateHapCode(count%2?this->hapcode:this->reduced,count%2?this->reduced:this->hapcode);
		}
		//std::cout<<"PhasedData:\n";
		//printArray(this->PhasedData->mGenotype);
		if (1 == numOfBlock ){
			this->hapcode=count%2?this->reduced:this->hapcode;
			break;
		}
		this->DataToPhase.reset(new SHEsisData(this->DataToPhase->getSampleNum(),numOfBlock,this->DataToPhase->getNumOfChrSet()));
		this->DataToPhase->vLabel=this->data->vLabel;
		this->DataToPhase->vQuantitativeTrait=this->data->vQuantitativeTrait;
		this->DataToPhase->mGenotype.resize(extension(this->PhasedData->mGenotype));
		this->DataToPhase->mGenotype=this->PhasedData->mGenotype;
		count++;
	};
	this->hap->PhasedData.resize(extension(this->PhasedData->mGenotype));
	this->hap->PhasedData=this->PhasedData->mGenotype;
	this->hap->SnpIdx.clear();
	this->hap->SnpIdx.push_back(0);
	this->hap->getResults();
	for(int i=0;i<this->hap->Results.haplotypes.size();i++){
		std::string ss;
		boost::shared_ptr<short[]> cur=this->hap->Results.haplotypes[i];
		boost::shared_ptr<short[]> res(new short[this->SnpIdx.size()]);
		int idx=0;
		for(int j=0;j<1;j++){
			hapMap::map_by<code>::const_iterator iter=this->hapcode.by<code>().find(cur[j]);
			BOOST_ASSERT(iter!=this->hapcode.by<code>().end());
			ss=iter->get<haplo>();
			std::cout<<"hap code:"<<cur[0]<<",haplotype:"<<ss<<"\n";
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
  		  	    	std::cerr<</*this->data->getallele*/(this->haplotypes[this->hap->Results.genotypes[i][j]][k]);
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
	//haplotypeBase:310
	this->hap->Results.haplotypes=this->haplotypes;
	this->hap->SnpIdx=this->SnpIdx;
	this->hap->data=this->data;
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
void HaplotypeLD::getHaplotypeSub(int start, int end,boost::shared_ptr<SHEsisData> d){
	std::vector<short> m(d->getSnpNum(),0);
	for(int i=start;i<=end;i++){
		m[this->SnpIdx[i]]=1;
	};
	this->hap.reset(new HaplotypeEM(d,end-start+1,m));
	this->hap->ShowResults(false);
	this->hap->startHaplotypeAnalysis();
}

void HaplotypeLD::updateData(boost::shared_ptr<SHEsisData> ret, boost::multi_array<short,3> phased) {
	ret->mGenotype=phased;
}

HaplotypeLD::~HaplotypeLD() {
	// TODO Auto-generated destructor stub
}

} /* namespace SHEsis */
