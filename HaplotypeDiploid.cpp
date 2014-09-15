/*
 * HaplotypeDiploid.cpp
 *
 *  Created on: Sep 9, 2014
 *      Author: ada
 */

#include "HaplotypeDiploid.h"

namespace SHEsis {

HaplotypeDiploid::HaplotypeDiploid(boost::shared_ptr<SHEsisData> data):HaplotypeBase(data),
		phased(1),PhasedData(boost::extents[data->getSampleNum()][data->getSnpNum()][2]),
		InterMediate(boost::extents[data->getSampleNum()][2][2]),err(0.00001){
BOOST_ASSERT(data->getNumOfChrSet() == 2);
for(int i=0;i<data->getSampleNum();i++){
	for(int j=0;j<2;j++){
		PhasedData[i][0][j]=data->mGenotype[i][0][j];
	};
};
}

HaplotypeDiploid::HaplotypeDiploid(boost::shared_ptr<SHEsisData>  data, int Snp, std::vector<short> mask):HaplotypeBase(data,mask),
		phased(1),PhasedData(boost::extents[data->getSampleNum()][data->getSnpNum()][2]),
		InterMediate(boost::extents[data->getSampleNum()][2][2]),err(0.00001){
BOOST_ASSERT(data->getNumOfChrSet() == 2);
for(int i=0;i<data->getSampleNum();i++){
	for(int j=0;j<2;j++){
		PhasedData[i][0][j]=data->mGenotype[i][0][j];
	};
};
}

HaplotypeDiploid::~HaplotypeDiploid() {
}

bool isMissing(boost::multi_array<short,3> data, int sample, int start, int end){
	BOOST_ASSERT(sample<data.shape()[0]);
	BOOST_ASSERT(start<=end);
	BOOST_ASSERT(end<data.shape()[1]);
	for(int i=start;i<=end;i++){
		for(int j=0;j<data.shape()[2];j++){
			if(0 == data[sample][i][j])
				return true;
		}
	}
	return false;
}

bool GenotypeEqual(boost::multi_array<short,3> data, int sample1, int sample2, int snp){
	if(((data[sample1][snp][0] == data[sample2][snp][0]) && (data[sample1][snp][1] == data[sample2][snp][1]))||
		((data[sample1][snp][1] == data[sample2][snp][0]) && (data[sample1][snp][0] == data[sample2][snp][1]))){
		return true;
	}else
		return false;
}

int getPhenotypeCode(std::vector<std::string>& v, std::string str){
	for(int i=0;i<v.size();i++){
		if(strcmp(v[i].c_str(),str.c_str()) == 0)
			return i;
	}
	v.push_back(str);
	return (v.size()-1);
}
void HaplotypeDiploid::ReturnGenotypeCode(int sample,short& geno1, short& geno2){
	std::stringstream p1("");
	std::stringstream p2("");
	geno1=-1;
	geno2=-1;
	for(int i=0;i<this->phased;i++){
		if(0 == this->PhasedData[sample][i][0]){
			geno1=0;
			break;
		}
		p1<<this->PhasedData[sample][i][0];
	};
	if(0 != geno1)
		geno1=getPhenotypeCode(this->InterMediateGenoCode,p1.str());

	for(int i=0;i<this->phased;i++){
		if(0 == this->PhasedData[sample][i][1]){
			geno2=0;
			break;
		}
		p2<<this->PhasedData[sample][i][1];
	};
	if(0 != geno2)
		geno2=getPhenotypeCode(this->InterMediateGenoCode,p2.str());
}

void HaplotypeDiploid::GenerateInterMediate(){
	this->InterMediateGenoCode.clear();
	InterMediateGenoCode.push_back("place holder");
	for(int i=0;i<data->getSampleNum();i++){
		short g1,g2;
		this->ReturnGenotypeCode(i,g1,g2);
		this->InterMediate[i][0][0]=g1;
		this->InterMediate[i][0][1]=g2;
		this->InterMediate[i][1][0]=this->data->mGenotype[i][this->phased][0];
		this->InterMediate[i][1][1]=this->data->mGenotype[i][this->phased][1];
	}
}


void HaplotypeDiploid::GenerateUniqueGenotype(){
	int idx=0;
	while(isMissing(this->InterMediate,idx, 0, 1))
		idx++;

	this->UniqueGenotypeIdx.push_back(idx++);
	this->UniqueGenotypeCount.push_back(1);
	for(;idx<this->data->getSampleNum();idx++){
		if(isMissing(this->InterMediate,idx,0,1))
			continue;
		for(int i=0;i<this->UniqueGenotypeIdx.size();i++){
			if(GenotypeEqual(this->InterMediate, idx,this->UniqueGenotypeIdx[i],0) &&
					GenotypeEqual(this->InterMediate, idx,this->UniqueGenotypeIdx[i],1)){
				this->UniqueGenotypeCount[i]++;
				break;
			}
			if(this->UniqueGenotypeIdx.size()-1 == i){
				this->UniqueGenotypeIdx.push_back(idx);
				this->UniqueGenotypeCount.push_back(1);
			}
		}
	}
	BOOST_ASSERT(this->UniqueGenotypeCount.size() == this->UniqueGenotypeIdx.size());
}

std::vector<short> getAlleleType(boost::multi_array<short,3> data, int sample,int snp){
	std::vector<short> alleleType;
	boost::unordered_map<short,int> hm;
	for(int i=0;i<2;i++){
		if(hm.end()==hm.find(data[sample][snp][i])){
			hm[data[sample][snp][i]]=1;
			alleleType.push_back(data[sample][snp][i]);
		}
	}
	return alleleType;
}

std::vector<short> getAlleleType(boost::multi_array<short,3> data,std::vector<int> sampleIdx,int snp){
	std::vector<short> alleleType;
	boost::unordered_map<short,int> hm;
	for(int k=0;k<sampleIdx.size();k++){
		for(int i=0;i<2;i++){
			if(hm.end()==hm.find(data[sampleIdx[k]][snp][i])){
				hm[data[sampleIdx[k]][snp][i]]=1;
				alleleType.push_back(data[sampleIdx[k]][snp][i]);
			}
		}
	}
	return alleleType;
}


int getHaploIdx(boost::shared_ptr<short[]> h, std::vector< boost::shared_ptr<short[]> > hap){
	for(int i=0;i<hap.size();i++){
		if(h[0]==hap[i][0] && h[1] == hap[i][1])
			return i;
	}
	return -1;
}

int PickTheOtherHaplo(boost::multi_array<short,3> data, int sample,std::vector< boost::shared_ptr<short[]> > hap,int idx1){
	boost::shared_ptr<short[]> p1=hap[idx1];
	boost::shared_ptr<short[]> p2=(new short[2]);
	if(p1[0] == data[sample][0][0]){
		p2[0]=data[sample][0][1];
	}else if(p1[0] == data[sample][0][1]){
		p2[0]=data[sample][0][0];
	}else{
		BOOST_ASSERT(0==1);
	}

	if(p1[1] == data[sample][1][0]){
		p2[1]=data[sample][1][1];
	}else if(p1[1] == data[sample][1][1]){
		p2[1]=data[sample][1][0];
	}else{
		BOOST_ASSERT(0==1);
	}
	return getHaploIdx(p2,hap);
}

bool ExistsHaploPair(HaploPair& ap, std::vector<HaploPair>& expanded){
	for(int i=0;i<expanded.size();i++){
		if(ap == expanded[i])
			return true;
	}
	return false;
}

void HaplotypeDiploid::generateAllPossibleHap(){
	std::vector<short> alleleType1=getAlleleType(this->InterMediate,this->UniqueGenotypeIdx,0);
	std::vector<short> alleleType2=getAlleleType(this->InterMediate,this->UniqueGenotypeIdx,1);
	boost::shared_ptr<short[]> hap;
	for(int i=0;i<alleleType1.size();i++){
		for(int j=0;j<alleleType2.size();j++){
			hap.reset(new short[2]);
			hap[0]=alleleType1[i];
			hap[1]=alleleType2[j];
			OneGenotypeExpandedHaplo::haploType.push_back(hap);
		}
	}
}
bool compatitable(boost::multi_array<short,3> data,int sample,boost::shared_ptr<short[]> hap){
	bool site1=false;
	bool site2=false;
	if(data[sample][0][0]==hap[0] || data[sample][0][0]==hap[1] || data[sample][0][1]==hap[0] || data[sample][0][1]==hap[1])
		site1=true;
	if(data[sample][1][0]==hap[0] || data[sample][1][0]==hap[1] || data[sample][1][1]==hap[0] || data[sample][1][1]==hap[1])
		site2=true;
	return (site1&&site2);
}

OneGenotypeExpandedHaplo HaplotypeDiploid::OneGenoExpandHaplo(int sample){
	OneGenotypeExpandedHaplo OneGenoHp;
	for(int i=0;i<OneGenoHp.haploType.size();i++){
		if(compatitable(this->InterMediate,sample,OneGenotypeExpandedHaplo::haploType[i])){
			HaploPair hp;
			hp.hap1=i;
			hp.hap2=PickTheOtherHaplo(this->InterMediate,sample,OneGenotypeExpandedHaplo::haploType,i);
			BOOST_ASSERT(-1 != hp.hap2);
			if(!ExistsHaploPair(hp,OneGenoHp.hp))
				OneGenoHp.hp.push_back(hp);
		};
	}
	return OneGenoHp;
}

void HaplotypeDiploid::ExpandAllGenotype(){
	this->Expanded.clear();
	for(int i=0;i<this->UniqueGenotypeIdx.size();i++){
		this->Expanded.push_back(
				OneGenoExpandHaplo(this->UniqueGenotypeIdx[i])
				);
	}
}
bool containsHap(HaploPair s,int hapidx){
	return (s.hap1 == hapidx || s.hap2 == hapidx);
}

void HaplotypeDiploid::CalculateFreq(){
	int hapcount=OneGenotypeExpandedHaplo::haploType.size();
	for(int i=0;i<hapcount;i++){
		OneGenotypeExpandedHaplo::hapfreq.push_back(1.0/(double)hapcount);
	};

	std::vector<double> H(OneGenotypeExpandedHaplo::haploType.size(),1.0/(double)hapcount);

	double E=1;
	while(sqrt(E)>this->err){
		std::vector<double> M(this->UniqueGenotypeCount.size(),0);
		for(int i=0;i<this->UniqueGenotypeCount.size();i++){
			for(int j=0;j<this->Expanded[i].hp.size();j++){
				int h1=this->Expanded[i].hp[j].hap1;
				int h2=this->Expanded[i].hp[j].hap2;
				M[i]+=H[h1]*H[h2];
			}
		}
		for(int n=0;n<hapcount;n++){
			double G=0;
			for(int j=0;j<this->UniqueGenotypeCount.size();j++){
				double U=0;
				for(int k=0;k<this->Expanded[j].hp.size();k++){
					if(containsHap(this->Expanded[j].hp[k],n)){
						U+=H[this->Expanded[j].hp[k].hap1]*H[this->Expanded[j].hp[k].hap2];
					}
				}
				G+=(double)this->UniqueGenotypeCount[j]/(double)this->data->getSampleNum()/2.0*(U/M[j]);
			}
		OneGenotypeExpandedHaplo::hapfreq[n]=G;
		}
		E=0;
		for(int n=0;n<hapcount;n++){
			E+=(OneGenotypeExpandedHaplo::hapfreq[n]-H[n])*(OneGenotypeExpandedHaplo::hapfreq[n]-H[n]);
			H[n]=OneGenotypeExpandedHaplo::hapfreq[n];
		}
	}

	for(int i=0;i<hapcount;i++){
		OneGenotypeExpandedHaplo::hapfreq[i]=H[i];
	}
}





} /* namespace SHEsis */
