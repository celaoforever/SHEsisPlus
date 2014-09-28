/*
 * HaplotypeDiploid.cpp
 *
 *  Created on: Sep 9, 2014
 *      Author: ada
 */

#include "HaplotypeEM.h"
#include "fisher.h"
#include "utility.h"

template
inline bool next_combination(std::vector<int>::iterator n_begin, std::vector<int>::iterator  n_end,
		std::vector<int>::iterator  r_begin, std::vector<int>::iterator  r_end);

namespace SHEsis {
std::vector<boost::shared_ptr<short[]> > OneGenotypeExpandedHaplo::haploType;
std::vector<double> OneGenotypeExpandedHaplo::hapfreq;
HaplotypeEM::HaplotypeEM(boost::shared_ptr<SHEsisData> data)
    : HaplotypeBase(data),
      phased(1),
      PhasedData(boost::extents[data->getSampleNum()][data->getSnpNum()][data->getNumOfChrSet()]),
      InterMediate(boost::extents[data->getSampleNum()][2][data->getNumOfChrSet()]),
      err(0.0001),
      missing(new bool[data->getSampleNum()]) {
  for (int i = 0; i < this->data->getSnpNum(); i++) {
    this->SnpIdx.push_back(i);
  }
  for (int i = 0; i < data->getSampleNum(); i++) {
    for (int j = 0; j < data->getNumOfChrSet(); j++) {
      PhasedData[i][0][j] = data->mGenotype[i][this->SnpIdx[0]][j];
    };
  };
  this->statMissing();
}

HaplotypeEM::HaplotypeEM(boost::shared_ptr<SHEsisData> data, int Snp,
                                   std::vector<short> mask)
    : HaplotypeBase(data, mask),
      phased(1),
      PhasedData(boost::extents[data->getSampleNum()][data->getSnpNum()][data->getNumOfChrSet()]),
      InterMediate(boost::extents[data->getSampleNum()][2][data->getNumOfChrSet()]),
      err(0.00001),
      missing(new bool[data->getSampleNum()]) {
  for (int i = 0; i < this->mask.size(); i++) {
    if (this->mask[i]) {
      this->SnpIdx.push_back(i);
    };
  };
  BOOST_ASSERT(this->SnpIdx.size() == Snp);

  for (int i = 0; i < data->getSampleNum(); i++) {
    for (int j = 0; j < data->getNumOfChrSet(); j++) {
      PhasedData[i][0][j] = data->mGenotype[i][this->SnpIdx[0]][j];
    };
  };
  this->statMissing();
}

HaplotypeEM::~HaplotypeEM() {}

void HaplotypeEM::sortGenotype(){
	for(int sample=0;sample<this->data->getSampleNum();sample++){
		for(int snp=0;snp<this->data->getSnpNum();snp++){
			std::vector tmp;
			for(int p=0;p<this->data->getNumOfChrSet();p++){
				tmp.push_back(this->data->mGenotype[sample][snp][p]);
			};
			std::sort(tmp.begin(),tmp.end());
			for(int p=0;p<this->data->getNumOfChrSet();p++){
				this->data->mGenotype[sample][snp][p]=tmp[p];
			}
		}
	}
}

void HaplotypeEM::statMissing() {
  for (int sample = 0; sample < this->data->getSampleNum(); sample++) {
    this->missing[sample] = false;
    for (int snp = 0; snp < this->SnpIdx.size(); snp++) {
    	for(int p=0;p<data->getNumOfChrSet();p++){
    		if(0 == this->data->mGenotype[sample][this->SnpIdx[snp]][p] ){
    			this->missing[sample] = true;
    			break;
    		}
    	}
		if(this->missing[sample])
			break;
    }
  }
}

bool isMissing(boost::multi_array<short, 3> data, int sample, int start,
               int end) {
  BOOST_ASSERT(sample < data.shape()[0]);
  BOOST_ASSERT(start <= end);
  BOOST_ASSERT(end < data.shape()[1]);
  for (int i = start; i <= end; i++) {
    for (int j = 0; j < data.shape()[2]; j++) {
      if (0 == data[sample][i][j]) return true;
    }
  }
  return false;
}

bool GenotypeEqual(boost::multi_array<short, 3> data, int sample1, int sample2,
                   int snp, int ploidy) {
	for(int p=0;p<ploidy;p++){
		if(data[sample1][snp][p]!=data[sample2][snp][p])
			return false;
	}
	return true;
}

int getPhenotypeCode(std::vector<std::string>& v, std::string str) {
  for (int i = 0; i < v.size(); i++) {
    if (strcmp(v[i].c_str(), str.c_str()) == 0) return i;
  }
  v.push_back(str);
  return (v.size() - 1);
}
void HaplotypeEM::ReturnGenotypeCode(int sample, std::vector<short>& geno) {
	if(missing[sample]){
		for(int p=0;p<this->data->getNumOfChrSet();p++){
			geno.push_back(0);
		}
		return;

	};
	for(int p=0;p>this->data->getNumOfChrSet();p++){
		std::stringstream p;
		for(int snp=0;snp<this->phased;snp++){
			p<<this->PhasedData[sample][snp][p];
		}
		geno.push_back(getPhenotypeCode(this->InterMediateGenoCode,p.str()));
	}
}

void HaplotypeEM::GenerateInterMediate() {
  this->InterMediateGenoCode.clear();
  InterMediateGenoCode.push_back("place holder");
  for (int sample = 0; sample < data->getSampleNum(); sample++) {
	  std::vector<short> geno;
	  this->ReturnGenotypeCode(sample,geno);
	  for(int p=0;p<this->data->getNumOfChrSet();p++){
		  this->InterMediate[sample][0][p]=geno[p];
		  this->InterMediate[sample][1][p]=this->data->mGenotype[sample][this->SnpIdx[this->phased]][p];
	  }
//    short g1, g2;
//    if (missing[i]) {
//      g1 = g2 = 0;
//    } else {
//      this->ReturnGenotypeCode(i, g1, g2);
//    }
//    this->InterMediate[i][0][0] = g1;
//    this->InterMediate[i][0][1] = g2;
//    this->InterMediate[i][1][0] =
//        this->data->mGenotype[i][this->SnpIdx[this->phased]][0];
//    this->InterMediate[i][1][1] =
//        this->data->mGenotype[i][this->SnpIdx[this->phased]][1];
  }
}

void HaplotypeEM::GenerateUniqueGenotype() {
  this->UniqueGenotypeIdx.clear();
  this->UniqueGenotypeCount.clear();
  this->Sample2Genotype.reset(new int[this->data->getSampleNum()]);
  int idx = 0;
  while (this->missing[idx]) {
    this->Sample2Genotype[idx] = -1;
    idx++;
  }

  this->UniqueGenotypeIdx.push_back(idx++);
  this->UniqueGenotypeCount.push_back(1);
  this->Sample2Genotype[idx - 1] = 0;
  for (; idx < this->data->getSampleNum(); idx++) {
    if (this->missing[idx]) {
      this->Sample2Genotype[idx] = -1;
      continue;
    }
    for (int i = 0; i < this->UniqueGenotypeIdx.size(); i++) {
      if (GenotypeEqual(this->InterMediate, idx, this->UniqueGenotypeIdx[i],
                        0, this->data->getNumOfChrSet()) &&
          GenotypeEqual(this->InterMediate, idx, this->UniqueGenotypeIdx[i],
                        1,this->data->getNumOfChrSet())) {
        this->UniqueGenotypeCount[i]++;
        this->Sample2Genotype[idx] = i;
        break;
      }
      if (this->UniqueGenotypeIdx.size() - 1 == i) {
        this->UniqueGenotypeIdx.push_back(idx);
        this->UniqueGenotypeCount.push_back(1);
        this->Sample2Genotype[idx] = this->UniqueGenotypeIdx.size() - 1;
        break;
      }
    }
  }
  BOOST_ASSERT(this->UniqueGenotypeCount.size() ==
               this->UniqueGenotypeIdx.size());
}

std::vector<short> getAlleleType(boost::multi_array<short, 3> data, int sample,
                                 int snp) {
  std::vector<short> alleleType;
  boost::unordered_map<short, int> hm;
  for (int i = 0; i < 2; i++) {
    if (hm.end() == hm.find(data[sample][snp][i])) {
      hm[data[sample][snp][i]] = 1;
      alleleType.push_back(data[sample][snp][i]);
    }
  }
  return alleleType;
}

std::vector<short> getAlleleType(boost::multi_array<short, 3> data,
                                 std::vector<int> sampleIdx, int snp,int ploidy) {
  std::vector<short> alleleType;
  boost::unordered_map<short, int> hm;
  for (int k = 0; k < sampleIdx.size(); k++) {
    for (int i = 0; i < ploidy; i++) {
      if (hm.end() == hm.find(data[sampleIdx[k]][snp][i])) {
        hm[data[sampleIdx[k]][snp][i]] = 1;
        alleleType.push_back(data[sampleIdx[k]][snp][i]);
      }
    }
  }
  return alleleType;
}

int getHaploIdx(boost::shared_ptr<short[]> h,
                std::vector<boost::shared_ptr<short[]> > hap) {
  for (int i = 0; i < hap.size(); i++) {
    if (h[0] == hap[i][0] && h[1] == hap[i][1]) return i;
  }
  return -1;
}

int PickTheOtherHaplo(boost::multi_array<short, 3> data, int sample,
                      std::vector<boost::shared_ptr<short[]> > hap, int idx1) {
  boost::shared_ptr<short[]> p1 = hap[idx1];
  boost::shared_ptr<short[]> p2(new short[2]);
  if (p1[0] == data[sample][0][0]) {
    p2[0] = data[sample][0][1];
  } else if (p1[0] == data[sample][0][1]) {
    p2[0] = data[sample][0][0];
  } else {
    BOOST_ASSERT(0 == 1);
  }

  if (p1[1] == data[sample][1][0]) {
    p2[1] = data[sample][1][1];
  } else if (p1[1] == data[sample][1][1]) {
    p2[1] = data[sample][1][0];
  } else {
    BOOST_ASSERT(0 == 1);
  }
  return getHaploIdx(p2, hap);
}

bool ExistsHaploCombination(HaploCombination& ap, std::vector<HaploCombination>& expanded) {
  for (int i = 0; i < expanded.size(); i++) {
    if (ap == expanded[i]) return true;
  }
  return false;
}

void HaplotypeEM::generateAllPossibleHap() {
  OneGenotypeExpandedHaplo::haploType.clear();
  std::vector<short> alleleType1 =
      getAlleleType(this->InterMediate, this->UniqueGenotypeIdx, 0);
  std::vector<short> alleleType2 =
      getAlleleType(this->InterMediate, this->UniqueGenotypeIdx, 1);
  boost::shared_ptr<short[]> hap;
  for (int i = 0; i < alleleType1.size(); i++) {
    for (int j = 0; j < alleleType2.size(); j++) {
      hap.reset(new short[2]);
      hap[0] = alleleType1[i];
      hap[1] = alleleType2[j];
      OneGenotypeExpandedHaplo::haploType.push_back(hap);
    }
  }
}

void HaplotypeEM::getCombination(){
	int ploidy=this->data->getNumOfChrSet();
	int total=ploidy*ploidy;
	//start from zero
	std::vector<int> c;
	std::vector<int> r;
	for(int i=0;i<ploidy;i++)
		c.push_back(i);

	for(int i=0;i<ploidy;i++)
		r.push_back(i);

	do{
		std::vector<int> tmp=r;
		this->combination.push_back(tmp);
	}while(next_combination(c.begin(),c.end(),r.begin(),r.end()));

}


bool GenotypeEqual(std::vector<short> vsnp, boost::multi_array<short, 3> originData,int sample,
                   int snp, int ploidy) {
	for(int p=0;p<ploidy;p++){
		if(vsnp[p]!=originData[sample][snp][p])
			return false;
	}
	return true;
}

OneGenotypeExpandedHaplo HaplotypeEM::oneGenoGetCombination(int sample){
	OneGenotypeExpandedHaplo res;
	std::vector< boost::shared_ptr<short[]> > v;
	for(int i=0;i<this->data->getNumOfChrSet();i++){
		for(int j=0;j<this->data->getNumOfChrSet();j++){
			boost::shared_ptr<short[]> hp(new short[2]);
			hp[0]=this->InterMediate[sample][0][i];
			hp[1]=this->InterMediate[sample][1][j];
			v.push_back(hp);
		}
	}

	for(int i=0;i<this->combination.size();i++){
		std::vector<int> cur=this->combination[i];
		std::vector<short> snp1;
		std::vector<short> snp2;
		for(int k=0;k<cur.size();k++){
			int curIdx=cur[k];
			snp1.push_back(v[curIdx][0]);
			snp2.push_back(v[curIdx][1]);
		}
		std::sort(snp1.begin(),snp1.end());
		std::sort(snp2.begin(),snp2.end());
		if(GenotypeEqual(snp1,this->InterMediate,sample,0,this->data->getNumOfChrSet()) &&
		   GenotypeEqual(snp2,this->InterMediate,sample,1,this->data->getNumOfChrSet())){
			HaploCombination hc;
			for(int n=0;n<cur.size();n++){
				hc.hapIdx.push_back(getHaploIdx(v[n],OneGenotypeExpandedHaplo::haploType));
			}
			if(!ExistsHaploCombination(hc,res.hp))
				res.hp.push_back(hc);
		}
	}
	return res;
}

bool compatitable(boost::multi_array<short, 3> data, int sample,
                  boost::shared_ptr<short[]> hap) {
  bool site1 = false;
  bool site2 = false;
  if (data[sample][0][0] == hap[0] || data[sample][0][1] == hap[0])
    site1 = true;
  if (data[sample][1][0] == hap[1] || data[sample][1][0] == hap[1])
    site2 = true;
  return (site1 && site2);
}

//OneGenotypeExpandedHaplo HaplotypeEM::OneGenoExpandHaplo(int sample) {
//  OneGenotypeExpandedHaplo OneGenoHp;
//  for (int i = 0; i < OneGenoHp.haploType.size(); i++) {
//    if (compatitable(this->InterMediate, sample,
//                     OneGenotypeExpandedHaplo::haploType[i])) {
//      HaploPair hp;
//      hp.hap1 = i;
//      hp.hap2 = PickTheOtherHaplo(this->InterMediate, sample,
//                                  OneGenotypeExpandedHaplo::haploType, i);
//      BOOST_ASSERT(-1 != hp.hap2);
//      if (!ExistsHaploPair(hp, OneGenoHp.hp)) OneGenoHp.hp.push_back(hp);
//    };
//  }
//  return OneGenoHp;
//}

void HaplotypeEM::ExpandAllGenotype() {
  this->Expanded.clear();
  for (int i = 0; i < this->UniqueGenotypeIdx.size(); i++) {
    this->Expanded.push_back(oneGenoGetCombination(this->UniqueGenotypeIdx[i]));
  }
}

bool containsHap(HaploCombination s, int hapidx) {
	for(int i=0;i<s.hapIdx.size();i++){
		if(s.hapIdx[i]==hapidx)
			return true;
	}
	return false;
}

void HaplotypeEM::CalculateFreq() {
  int hapcount = OneGenotypeExpandedHaplo::haploType.size();
  for (int i = 0; i < hapcount; i++) {
    OneGenotypeExpandedHaplo::hapfreq.push_back(1.0 / (double)hapcount);
  };

  std::vector<double> H(OneGenotypeExpandedHaplo::haploType.size(),
                        1.0 / (double)hapcount);

  double E = 1;
  while (E > this->err) {
    //std::vector<double> M(this->UniqueGenotypeCount.size(), 0);
	  boost::shared_ptr<double[]> M(new double[this->UniqueGenotypeCount.size()]);
    for (int i = 0; i < this->UniqueGenotypeCount.size(); i++) {
      for (int j = 0; j < this->Expanded[i].hp.size(); j++) {
    	  double m=1;
    	  for(int k=0;k<this->data->getNumOfChrSet();k++){
    		  m*=this->Expanded[i].hp[j].hapIdx[k];
    	  }
        M[i] += m;
      }
    }
    for (int n = 0; n < hapcount; n++) {
      double G = 0;
      for (int j = 0; j < this->UniqueGenotypeCount.size(); j++) {
        double U = 0;
        for (int k = 0; k < this->Expanded[j].hp.size(); k++) {
          if (containsHap(this->Expanded[j].hp[k], n)) {
        	  double u=1;
        	  for(int kk=0;kk<this->data->getNumOfChrSet();kk++){
        		  u*=H[this->Expanded[j].hp[k].hapIdx[kk]];
        	  }
        	  U+=u;
          }
        }
        G += (double)this->UniqueGenotypeCount[j] /
             (double)this->data->getSampleNum() / (double)this->data->getNumOfChrSet() * (U / M[j]);
      }
      OneGenotypeExpandedHaplo::hapfreq[n] = G;
    }
    E = 0;
    for (int n = 0; n < hapcount; n++) {
      E += (OneGenotypeExpandedHaplo::hapfreq[n] - H[n]) *
           (OneGenotypeExpandedHaplo::hapfreq[n] - H[n]);
      H[n] = OneGenotypeExpandedHaplo::hapfreq[n];
    }
  }
}

void HaplotypeEM::getFinalHap() {
  for (int i = 0; i < this->UniqueGenotypeCount.size(); i++) {
    this->Expanded[i].freq = 0;
    for (int j = 0; j < this->Expanded[i].hp.size(); j++) {
    	double freq=1;
    	for(int k=0;k<this->data->getNumOfChrSet();k++){
    		int h=this->Expanded[i].hp[j].hapIdx[k];
    		freq*=OneGenotypeExpandedHaplo::hapfreq[h];
    	};
      if (this->Expanded[i].freq < freq) {
    	  this->Expanded[i].finalHapIdxInHp=j;
    	  this->Expanded[i].freq = freq;
      }
    }
  }
}

int locate(boost::multi_array<short,3> inter,int  sample, int ploidy,
		short allele,std::vector<bool>& isPhased){
	for(int p=0;p<ploidy;p++){
		if(inter[sample][0][p]==allele && !isPhased[p]){
			isPhased[p]=true;
			return p;
		}
	}
	return -1;
}

std::vector<int> HaplotypeEM::getSampleIdx(int genotype) {
  std::vector<int> res;
  for (int i = 0; i < this->data->getSampleNum(); i++) {
    if (this->Sample2Genotype[i] == genotype) res.push_back(i);
  }
  return res;
}
void HaplotypeEM::PhaseCurrent() {
  for (int i = 0; i < this->UniqueGenotypeCount.size(); i++) {
    std::vector<int> samples = getSampleIdx(i);

    for(int sampleIdx=0;sampleIdx<samples.size();sampleIdx++){
    	int finalHapIdx=this->Expanded[i].finalHapIdxInHp;
    	HaploCombination CurCombination=this->Expanded[i].hp[finalHapIdx];
    	std::vector<bool> isPhased;
    	isPhased.resize(this->data->getNumOfChrSet(),false);
    	for(int p=0;p<this->data->getNumOfChrSet();p++){
    		short CurHapAllele1=OneGenotypeExpandedHaplo::haploType[CurCombination[p]][0];
    		short CurHapAllele2=OneGenotypeExpandedHaplo::haploType[CurCombination[p]][1];
    		int whichPloidy=locate(this->InterMediate,sampleIdx,this->data->getNumOfChrSet(),CurHapAllele1,isPhased);
    		BOOST_ASSERT(whichPloidy!=-1);
    		this->PhasedData[sampleIdx][this->phased][whichPloidy]=CurHapAllele2;
    	}

    }
//    int h1 = this->Expanded[i].finalhap1;
//    int h2 = this->Expanded[i].finalhap2;
//    short p00 = OneGenotypeExpandedHaplo::haploType[h1][0];
//    short p01 = OneGenotypeExpandedHaplo::haploType[h1][1];
//    short p10 = OneGenotypeExpandedHaplo::haploType[h2][0];
//    short p11 = OneGenotypeExpandedHaplo::haploType[h2][1];
//    for (int j = 0; j < samples.size(); j++) {
//      int idx = samples[j];
//      short q00 = this->InterMediate[idx][0][0];
//      short q01 = this->InterMediate[idx][1][0];
//      short q10 = this->InterMediate[idx][0][1];
//      short q11 = this->InterMediate[idx][1][1];
//      if (p00 == q00 && p10 == q10) {
//        BOOST_ASSERT((p01 == q01 && p11 == q11) || (p01 == q11 && p11 == q01));
//        this->PhasedData[idx][this->phased][0] = p01;
//        this->PhasedData[idx][this->phased][1] = p11;
//      } else if (p00 == q10 && p10 == q00) {
//        BOOST_ASSERT((p01 == q01 && p11 == q11) || (p01 == q11 && p11 == q01));
//        this->PhasedData[idx][this->phased][0] = p11;
//        this->PhasedData[idx][this->phased][1] = p01;
//
//      } else {
//        BOOST_ASSERT(1 == 0);
//      }
//    }
  }
  this->phased++;
};

void HaplotypeEM::getResults() {
  boost::unordered_map<std::string, int> hm;
  boost::shared_ptr<short[]> haplo;
  int idx = 0;
  std::stringstream p1;
  std::stringstream p2;
  for (int sample = 0; sample < this->data->getSampleNum(); sample++) {
    if (this->missing[sample]) {
      this->Results.genotypes[sample][0] = -1;
      this->Results.genotypes[sample][1] = -1;
      continue;
    }
    p1.str("");
    p2.str("");
    for (int snp = 0; snp < this->SnpIdx.size(); snp++) {
      p1 << this->PhasedData[sample][snp][0] << ",";
      p2 << this->PhasedData[sample][snp][1] << ",";
    }
    if (hm.find(p1.str()) == hm.end()) {
      haplo.reset(new short[this->SnpIdx.size()]);
      for (int snp = 0; snp < this->SnpIdx.size(); snp++) {
        haplo[snp] = this->PhasedData[sample][snp][0];
      }
      this->Results.haplotypes.push_back(haplo);
      hm[p1.str()] = idx++;
      BOOST_ASSERT(this->Results.haplotypes.size() == idx);
    };
    if (hm.find(p2.str()) == hm.end()) {
      haplo.reset(new short[this->SnpIdx.size()]);
      for (int snp = 0; snp < this->SnpIdx.size(); snp++) {
        haplo[snp] = this->PhasedData[sample][snp][1];
      }
      this->Results.haplotypes.push_back(haplo);

      hm[p2.str()] = idx++;
      BOOST_ASSERT(this->Results.haplotypes.size() == idx);
    };
    this->Results.genotypes[sample][0] = hm[p1.str()];
    this->Results.genotypes[sample][1] = hm[p2.str()];
  };

  int haploNum = hm.size();
  this->Results.CaseCount.reset(new int[haploNum]);
  this->Results.ControlCount.reset(new int[haploNum]);
  for (int i = 0; i < this->Results.haplotypes.size(); i++) {
    this->Results.CaseCount[i] = 0;
    this->Results.ControlCount[i] = 0;
  }

  for (int sample = 0; sample < this->data->getSampleNum(); sample++) {
    if (this->missing[sample]) continue;
    if (CASE == this->data->vLabel[sample]) {
      int h1 = this->Results.genotypes[sample][0];
      int h2 = this->Results.genotypes[sample][1];
      this->Results.CaseCount[h1]++;
      this->Results.CaseCount[h2]++;
    } else if (CONTROL == this->data->vLabel[sample]) {
      int h1 = this->Results.genotypes[sample][0];
      int h2 = this->Results.genotypes[sample][1];
      this->Results.ControlCount[h1]++;
      this->Results.ControlCount[h2]++;
    }
  };
}

void HaplotypeEM::startHaplotypeAnalysis() {
  while (this->phased < this->SnpIdx.size()) {
    if (!this->silent && this->phased % 10 == 0) {
      int per = 100 * (double)this->phased / (double)this->SnpIdx.size();
      printf("\rProgress:%d%%", per);
      fflush(stdout);
    }
    this->GenerateInterMediate();
    this->GenerateUniqueGenotype();
    this->generateAllPossibleHap();
    this->ExpandAllGenotype();
    this->CalculateFreq();
    this->getFinalHap();
    this->PhaseCurrent();
  };
  if (!this->silent) {
    printf("\rProgress:%d%%\n", 100);
    fflush(stdout);
  }

  this->getResults();
}

} /* namespace SHEsis */
