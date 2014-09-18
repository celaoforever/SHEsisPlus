/*
 * SHEsisData.cpp
 *
 *  Created on: Aug 3, 2014
 *      Author: ionadmin
 */

#include "SHEsisData.h"
#include <algorithm>
#include <boost/assert.hpp>
#include <sstream>
#include <iostream>
#include <boost/foreach.hpp>
namespace SHEsis {


SHEsisData::SHEsisData(int SampleNum, int SnpNum, int NumOfChrSet)
:SampleNum(SampleNum),SnpNum(SnpNum),NumOfChrSet(NumOfChrSet),
 mGenotype(boost::extents[SampleNum][SnpNum][NumOfChrSet]),
 vLocusInfo(SnpNum),CaseNum(-1),ControlNum(-1),codeIdx(1),
 vLabel(SampleNum),vPermutateLabel(SampleNum)
{
};


SHEsisData::~SHEsisData() {
	this->vLabel.clear();
	this->vPermutateLabel.clear();
	this->vLocusInfo.clear();
}

std::string getStrFromSortedGenotype(std::vector<short> v){
	std::stringstream ss;
	for(int i=0;i<v.size()-1;i++)
	{
		ss<<v[i]<<"/";
	}
	ss<<v[v.size()-1];
	return ss.str();
}

short SHEsisData::GetAlleleCode(std::string const val){
	boost::unordered_map<short,std::string>::iterator iter;
	for(iter=this->code2allele.begin();iter!=this->code2allele.end();iter++){
		if(val == iter->second)
			return iter->first;
	}
	this->code2allele[this->codeIdx]=val;
	return this->codeIdx++;
}

void SHEsisData::printLocusInfo()
{
	for(int i=0;i<this->vLocusInfo.size();i++)
	{
		boost::unordered_map<short, double> ::iterator map_it;
		std::cout<<"\nLocus "<<i;
		std::cout<<"\nCaseAlleleCount:\n";
		for(map_it=this->vLocusInfo[i].CaseAlleleCount.begin();map_it!=this->vLocusInfo[i].CaseAlleleCount.end();map_it++){
			std::cout<<map_it->first<<":"<<map_it->second<<", ";
		}
		std::cout<<"\nControlAlleleCount:\n";
		for(map_it=this->vLocusInfo[i].ControlAlleleCount.begin();map_it!=this->vLocusInfo[i].ControlAlleleCount.end();map_it++){
			std::cout<<map_it->first<<":"<<map_it->second<<", ";
		}
		std::cout<<"\nBothAlleleCount:\n";
		for(map_it=this->vLocusInfo[i].BothAlleleCount.begin();map_it!=this->vLocusInfo[i].BothAlleleCount.end();map_it++){
			std::cout<<map_it->first<<":"<<map_it->second<<", ";
		}
		boost::unordered_map<std::string, double> ::iterator map_it2;
		std::cout<<"\nCaseGenotypeCount:\n";
		for(map_it2=this->vLocusInfo[i].CaseGenotypeCount.begin();map_it2!=this->vLocusInfo[i].CaseGenotypeCount.end();map_it2++){
			std::cout<<map_it2->first<<":"<<map_it2->second<<", ";
		}
		std::cout<<"\nControlGenotypeCount:\n";
		for(map_it2=this->vLocusInfo[i].ControlGenotypeCount.begin();map_it2!=this->vLocusInfo[i].ControlGenotypeCount.end();map_it2++){
			std::cout<<map_it2->first<<":"<<map_it2->second<<", ";
		}
		std::cout<<"\nBothGenotypeCount:\n";
		for(map_it2=this->vLocusInfo[i].BothGenotypeCount.begin();map_it2!=this->vLocusInfo[i].BothGenotypeCount.end();map_it2++){
			std::cout<<map_it2->first<<":"<<map_it2->second<<", ";
		}
	}
}

void SHEsisData::getCaseAndControlNum(){
	this->CaseNum=0;this->ControlNum=0;
	for(int iSample=0;iSample<this->SampleNum;iSample++){
		if(CASE == this->vLabel[iSample]){
			this->CaseNum++;
		}else if(CONTROL == this->vLabel[iSample])
		{
			this->ControlNum++;
		}
	}
}

int SHEsisData::getCaseNum(){
	if(-1 == this->CaseNum)
		this->getCaseAndControlNum();
	return this->CaseNum;
}

int SHEsisData::getControlNum(){
	if(-1 == this->ControlNum)
		this->getCaseAndControlNum();
	return this->ControlNum;
}

void SHEsisData::statCount(std::vector< SampleStatus > & label){
	std::vector<short> geno;
	for(int iSample=0;iSample<this->SampleNum;iSample++){
		for(int iSnp=0;iSnp<this->SnpNum;iSnp++){
			geno.clear();
			for(int iChrset=0;iChrset<this->NumOfChrSet;iChrset++){
				BOOST_ASSERT(iSample<this->mGenotype.shape()[0] &&
						iSnp<this->mGenotype.shape()[1] &&
						iChrset<this->mGenotype.shape()[2]);
				short Cur=this->mGenotype[iSample][iSnp][iChrset];
				if(GENOTYPE_MISSING ==Cur)
					continue;
				if(CASE == label[iSample] ){
					if(this->vLocusInfo[iSnp].CaseAlleleCount.end() == this->vLocusInfo[iSnp].CaseAlleleCount.find(Cur))
						this->vLocusInfo[iSnp].CaseAlleleCount[Cur]=1;
					else
						(this->vLocusInfo[iSnp].CaseAlleleCount[Cur])++;
					if(this->vLocusInfo[iSnp].ControlAlleleCount.end() == this->vLocusInfo[iSnp].ControlAlleleCount.find(Cur))
						this->vLocusInfo[iSnp].ControlAlleleCount[Cur]=0;
					if(iChrset == (this->NumOfChrSet-1)){
						geno.push_back(Cur);
						if(this->NumOfChrSet > geno.size())
							continue;   //indicating there is missing data within this site for the specific individual
						std::sort(geno.begin(),geno.end());
						std::string genoStr=getStrFromSortedGenotype(geno);
						if(this->vLocusInfo[iSnp].CaseGenotypeCount.end() == this->vLocusInfo[iSnp].CaseGenotypeCount.find(genoStr))
							this->vLocusInfo[iSnp].CaseGenotypeCount[genoStr]=1;
						else
							(this->vLocusInfo[iSnp].CaseGenotypeCount[genoStr])++;

						if(this->vLocusInfo[iSnp].ControlGenotypeCount.end() == this->vLocusInfo[iSnp].ControlGenotypeCount.find(genoStr))
							this->vLocusInfo[iSnp].ControlGenotypeCount[genoStr]=0;
					}else
						geno.push_back(Cur);

					continue;
				};
				if(CONTROL == label[iSample]){
					if(this->vLocusInfo[iSnp].ControlAlleleCount.end() == this->vLocusInfo[iSnp].ControlAlleleCount.find(Cur))
						this->vLocusInfo[iSnp].ControlAlleleCount[Cur]=1;
					else
						(this->vLocusInfo[iSnp].ControlAlleleCount[Cur])++;
					if(this->vLocusInfo[iSnp].CaseAlleleCount.end() == this->vLocusInfo[iSnp].CaseAlleleCount.find(Cur))
						this->vLocusInfo[iSnp].CaseAlleleCount[Cur]=0;
					if(iChrset == (this->NumOfChrSet-1)){
						geno.push_back(Cur);
						if(this->NumOfChrSet > geno.size())
							continue;   //indicating there is missing data within this site for the specific individual
						std::sort(geno.begin(),geno.end());
						std::string genoStr=getStrFromSortedGenotype(geno);
						if(this->vLocusInfo[iSnp].ControlGenotypeCount.end() == this->vLocusInfo[iSnp].ControlGenotypeCount.find(genoStr))
							this->vLocusInfo[iSnp].ControlGenotypeCount[genoStr]=1;
						else
							(this->vLocusInfo[iSnp].ControlGenotypeCount[genoStr])++;
						if(this->vLocusInfo[iSnp].CaseGenotypeCount.end() == this->vLocusInfo[iSnp].CaseGenotypeCount.find(genoStr))
							this->vLocusInfo[iSnp].CaseGenotypeCount[genoStr]=0;
					}else
						geno.push_back(Cur);
				}

			}

		}
	}
	for(int iSnp=0;iSnp<this->SnpNum;iSnp++){
		BOOST_ASSERT(this->vLocusInfo[iSnp].CaseGenotypeCount.size()
				== this->vLocusInfo[iSnp].ControlGenotypeCount.size());
		boost::unordered_map<std::string, double>::iterator iter;
		for(iter=this->vLocusInfo[iSnp].CaseGenotypeCount.begin();
				iter!=this->vLocusInfo[iSnp].CaseGenotypeCount.end();iter++){
			this->vLocusInfo[iSnp].BothGenotypeCount[iter->first]=
					this->vLocusInfo[iSnp].CaseGenotypeCount[iter->first]+
					this->vLocusInfo[iSnp].ControlGenotypeCount[iter->first];
		}

		BOOST_ASSERT(this->vLocusInfo[iSnp].CaseAlleleCount.size()
				== this->vLocusInfo[iSnp].ControlAlleleCount.size());
		boost::unordered_map<short,double>::iterator iter2;
		for(iter2=this->vLocusInfo[iSnp].CaseAlleleCount.begin();
				iter2!=this->vLocusInfo[iSnp].CaseAlleleCount.end();iter2++){
			this->vLocusInfo[iSnp].BothAlleleCount[iter2->first]=
					this->vLocusInfo[iSnp].CaseAlleleCount[iter2->first]+
					this->vLocusInfo[iSnp].ControlAlleleCount[iter2->first];
		}
	}


}

}
