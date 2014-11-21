/*
 * HaplotypeDiploid.cpp
 *
 *  Created on: Sep 9, 2014
 *      Author: ada
 */

#include "HaplotypeEM.h"
#include "fisher.h"
#include "utility.h"

template bool next_combination(std::vector<int>::iterator n_begin,
                               std::vector<int>::iterator n_end,
                               std::vector<int>::iterator r_begin,
                               std::vector<int>::iterator r_end);

namespace SHEsis {
std::vector<boost::shared_ptr<short[]> > OneGenotypeExpandedHaplo::haploType;
std::vector<double> OneGenotypeExpandedHaplo::hapfreq;
HaplotypeEM::HaplotypeEM(boost::shared_ptr<SHEsisData> data)
    : HaplotypeBase(data),
      phased(1),
      PhasedData(boost::extents[data->getSampleNum()][data->getSnpNum()]
                               [data->getNumOfChrSet()]),
      InterMediate(
          boost::extents[data->getSampleNum()][2][data->getNumOfChrSet()]),
      err(0.0001),
      missing(new bool[data->getSampleNum()]) {
  //	this->sortGenotype();
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
      PhasedData(boost::extents[data->getSampleNum()][data->getSnpNum()]
                               [data->getNumOfChrSet()]),
      InterMediate(
          boost::extents[data->getSampleNum()][2][data->getNumOfChrSet()]),
      err(0.00001),
      missing(new bool[data->getSampleNum()]) {
  //	this->sortGenotype();
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

void HaplotypeEM::sortGenotype() {
  for (int sample = 0; sample < this->data->getSampleNum(); sample++) {
    for (int snp = 0; snp < this->data->getSnpNum(); snp++) {
      std::vector<short> tmp;
      for (int p = 0; p < this->data->getNumOfChrSet(); p++) {
        tmp.push_back(this->data->mGenotype[sample][snp][p]);
      };
      std::sort(tmp.begin(), tmp.end());
      for (int p = 0; p < this->data->getNumOfChrSet(); p++) {
        this->data->mGenotype[sample][snp][p] = tmp[p];
      }
    }
  }

  std::cout << "sorted:\n";
  for (int sample = 0; sample < this->data->getSampleNum(); sample++) {
    for (int snp = 0; snp < this->data->getSnpNum(); snp++) {
      for (int p = 0; p < this->data->getNumOfChrSet(); p++) {
        std::cout << this->data->mGenotype[sample][snp][p] << "/";
      }
      std::cout << " ";
    }
    std::cout << "\n";
  };
}

void HaplotypeEM::sortInterMediate() {
  for (int sample = 0; sample < this->data->getSampleNum(); sample++) {
    for (int snp = 0; snp < 2; snp++) {
      std::vector<short> tmp;
      for (int p = 0; p < this->data->getNumOfChrSet(); p++) {
        tmp.push_back(this->InterMediate[sample][snp][p]);
      };
      std::sort(tmp.begin(), tmp.end());
      for (int p = 0; p < this->data->getNumOfChrSet(); p++) {
        this->InterMediate[sample][snp][p] = tmp[p];
      }
    }
  }

  std::cout << "sorted:\n";
  for (int sample = 0; sample < this->data->getSampleNum(); sample++) {
    for (int snp = 0; snp < 2; snp++) {
      for (int p = 0; p < this->data->getNumOfChrSet(); p++) {
        std::cout << this->InterMediate[sample][snp][p] << "/";
      }
      std::cout << " ";
    }
    std::cout << "\n";
  };
}

void HaplotypeEM::statMissing() {
  for (int sample = 0; sample < this->data->getSampleNum(); sample++) {
    this->missing[sample] = false;
    for (int snp = 0; snp < this->SnpIdx.size(); snp++) {
      for (int p = 0; p < data->getNumOfChrSet(); p++) {
        if (0 == this->data->mGenotype[sample][this->SnpIdx[snp]][p]) {
          this->missing[sample] = true;
          if (CASE == this->data->vLabel[sample])
            this->NonmissingCase++;
          else if (CONTROL == this->data->vLabel[sample])
            this->NonmissingCtrl++;
          else if (this->data->vQuantitativeTrait.size())
            this->NonmissingCase++;
          break;
        }
      }
      if (this->missing[sample]) break;
    }
  }
  if (this->data->vQuantitativeTrait.size())
    this->NonmissingCase = this->data->getSampleNum() - this->NonmissingCase;
  else {
    this->NonmissingCase = this->data->getCaseNum() - this->NonmissingCase;
    this->NonmissingCtrl = this->data->getControlNum() - this->NonmissingCtrl;
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
  std::vector<short> vSnp1;
  std::vector<short> vSnp2;

  for (int p = 0; p < ploidy; p++) {
    vSnp1.push_back(data[sample1][snp][p]);
    vSnp2.push_back(data[sample2][snp][p]);
  }
  std::sort(vSnp1.begin(), vSnp1.end());
  std::sort(vSnp2.begin(), vSnp2.end());
  for (int p = 0; p < ploidy; p++) {
    if (vSnp1[p] != vSnp2[p]) return false;
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
  if (missing[sample]) {
    for (int p = 0; p < this->data->getNumOfChrSet(); p++) {
      geno.push_back(0);
    }
    return;
  };
  for (int p = 0; p < this->data->getNumOfChrSet(); p++) {
    std::stringstream pp;
    for (int snp = 0; snp < this->phased; snp++) {
      pp << this->PhasedData[sample][snp][p];
    }
    geno.push_back(getPhenotypeCode(this->InterMediateGenoCode, pp.str()));
  }
}

void HaplotypeEM::GenerateInterMediate() {
  this->InterMediateGenoCode.clear();
  InterMediateGenoCode.push_back("place holder");
  for (int sample = 0; sample < data->getSampleNum(); sample++) {
    std::vector<short> geno;
    this->ReturnGenotypeCode(sample, geno);
    for (int p = 0; p < this->data->getNumOfChrSet(); p++) {
      this->InterMediate[sample][0][p] = geno[p];
      this->InterMediate[sample][1][p] =
          this->data->mGenotype[sample][this->SnpIdx[this->phased]][p];
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
      if (GenotypeEqual(this->InterMediate, idx, this->UniqueGenotypeIdx[i], 0,
                        this->data->getNumOfChrSet()) &&
          GenotypeEqual(this->InterMediate, idx, this->UniqueGenotypeIdx[i], 1,
                        this->data->getNumOfChrSet())) {
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
  //  std::cout<<"UniqieGenotypeIdx:";
  //  for(int i=0;i<this->UniqueGenotypeIdx.size();i++){
  //	  std::cout<<this->UniqueGenotypeIdx[i]<<",";
  //  }
  //  std::cout<<"\n";
  //  std::cout<<"Sample2Genotype:";
  //  for(int i=0;i<this->data->getSampleNum();i++){
  //	  std::cout<<this->Sample2Genotype[i]<<",";
  //  }
  //  std::cout<<"\n";
  //  std::cout<<"UniqeGenotypeCount.Size:"<<this->UniqueGenotypeCount.size()<<"\n";
  BOOST_ASSERT(this->UniqueGenotypeCount.size() ==
               this->UniqueGenotypeIdx.size());
}

std::vector<short> getAlleleType(boost::multi_array<short, 3> data, int sample,
                                 int snp, int ploidy) {
  std::vector<short> alleleType;
  boost::unordered_map<short, int> hm;
  for (int i = 0; i < ploidy; i++) {
    if (hm.end() == hm.find(data[sample][snp][i])) {
      hm[data[sample][snp][i]] = 1;
      alleleType.push_back(data[sample][snp][i]);
    }
  }
  return alleleType;
}

std::vector<short> getAlleleType(boost::multi_array<short, 3> data,
                                 std::vector<int> sampleIdx, int snp,
                                 int ploidy) {
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
    bool equal = true;
    for (int p = 0; p < 2; p++) {
      equal = equal && (hap[i][p] == h[p]);
    }
    if (equal) return i;
  }
  return -1;
}

bool ExistsHaploCombination(HaploCombination& ap,
                            std::vector<HaploCombination>& expanded) {
  //	std::cout<<"expanded:";
  //	for(int i=0;i<expanded.size();i++){
  //		for(int j=0;j<expanded[i].hapIdx.size();j++){
  //			std::cout<<expanded[i].hapIdx[j]<<",";
  //		};
  //		std::cout<<" ";
  //	}
  //	std::cout<<"\ncurrent combination:";
  //	for(int k=0;k<ap.hapIdx.size();k++){
  //		std::cout<<ap.hapIdx[k]<<",";
  //	}
  //	std::cout<<"\n";
  std::sort(ap.hapIdx.begin(), ap.hapIdx.end());
  for (int i = 0; i < expanded.size(); i++) {
    if (ap == expanded[i]) return true;
  }
  return false;
}

void HaplotypeEM::generateAllPossibleHap() {
  OneGenotypeExpandedHaplo::haploType.clear();
  OneGenotypeExpandedHaplo::hapfreq.clear();
  std::vector<short> alleleType1 =
      getAlleleType(this->InterMediate, this->UniqueGenotypeIdx, 0,
                    this->data->getNumOfChrSet());
  std::vector<short> alleleType2 =
      getAlleleType(this->InterMediate, this->UniqueGenotypeIdx, 1,
                    this->data->getNumOfChrSet());
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

void HaplotypeEM::getCombination() {
  int ploidy = this->data->getNumOfChrSet();
  int total = ploidy * ploidy;
  // start from zero
  std::vector<int> c;
  std::vector<int> r;
  for (int i = 0; i < total; i++) c.push_back(i);

  for (int i = 0; i < ploidy; i++) r.push_back(i);

  do {
    std::vector<int> tmp = r;
    this->combination.push_back(tmp);
    //		std::cout<<"combination "<<this->combination.size()<<":";
    //		for(int i=0;i<tmp.size();i++){
    //			std::cout<<tmp[i]<<",";
    //		}
    //		std::cout<<"\n";
  } while (next_combination(c.begin(), c.end(), r.begin(), r.end()));
}

bool GenotypeEqual(std::vector<short> vsnp,
                   boost::multi_array<short, 3> originData, int sample, int snp,
                   int ploidy) {
  std::vector<short> vsnp1;
  for (int p = 0; p < ploidy; p++) {
    vsnp1.push_back(originData[sample][snp][p]);
  }
  std::sort(vsnp1.begin(), vsnp1.end());
  for (int p = 0; p < ploidy; p++) {
    if (vsnp[p] != vsnp1[p]) return false;
  }
  return true;
}

bool ExistsHap(boost::shared_ptr<short[]> hap,
               std::vector<boost::shared_ptr<short[]> >& v) {
  for (int i = 0; i < v.size(); i++) {
    if (hap[0] == v[i][0] && hap[1] == v[i][1]) return true;
  }
  return false;
}
bool containsRepeatHap(std::vector<int>& com, std::vector<int>& repeat) {
  std::cout << "\ncom:\n";
  for (int i = 0; i < com.size(); i++) {
    std::cout << com[i] << ",";
  }
  std::cout << "\nrepeat:\n";
  for (int i = 0; i < repeat.size(); i++) {
    std::cout << repeat[i] << ",";
  }
  std::cout << "\n";
  for (int i = 0; i < com.size(); i++) {
    for (int j = 0; j < repeat.size(); j++) {
      if (com[i] == repeat[j]) return true;
    }
  }
  return false;
}

long choose(boost::shared_ptr<int[]> got, int n_chosen, int len, int at,
            int max_types, std::vector<std::vector<int> >& combination) {
  int i;
  long count = 0;
  if (n_chosen == len) {
    if (!got) return 1;
    std::vector<int> curCom;
    for (i = 0; i < len; i++) {
      curCom.push_back(got[i]);
    }
    combination.push_back(curCom);
    // printf("%s\t", donuts[got[i]]);
    // printf("\n");
    return 1;
  }

  for (i = at; i < max_types; i++) {
    if (got) got[n_chosen] = i;
    count += choose(got, n_chosen + 1, len, i, max_types, combination);
  }
  return count;
}

void OneSampleGetCombination(
    std::vector<boost::shared_ptr<short[]> >& Haplotypes, int ploidy,
    std::vector<std::vector<int> >& combination) {
  boost::shared_ptr<int[]> got(new int[ploidy + 1]);
  choose(got, 0, ploidy, 0, Haplotypes.size(), combination);
  //	std::cout<<"Combination.size:"<<combination.size()<<"\n";
  //	std::cout<<"One sample haplotypes:";
  //	for(int i=0;i<Haplotypes.size();i++){
  //		std::cout<<Haplotypes[i][0]<<Haplotypes[i][1]<<",";
  //	}
  //	std::cout<<"\nCombinatin:";
  //	for(int i=0;i<combination.size();i++){
  //		for(int j=0;j<combination[i].size();j++){
  //			int hap=combination[i][j];
  //			std::cout<<Haplotypes[hap][0]<<Haplotypes[hap][1]<<",";
  //			//std::cout<<combination[i][j]<<",";
  //		}
  //		std::cout<<"\n";
  //	}
}

OneGenotypeExpandedHaplo HaplotypeEM::oneGenoGetCombination(int sample) {
  OneGenotypeExpandedHaplo res;
  std::vector<boost::shared_ptr<short[]> > OneSampleHaplotypes;
  boost::shared_ptr<short[]> curHap;
  std::vector<short> alleleType1 = getAlleleType(this->InterMediate, sample, 0,
                                                 this->data->getNumOfChrSet());
  std::vector<short> alleleType2 = getAlleleType(this->InterMediate, sample, 1,
                                                 this->data->getNumOfChrSet());
  //	std::cout<<"Current genotype:";
  //	for(int k=0;k<2;k++){
  //		for(int i=0;i<this->data->getNumOfChrSet();i++){
  //			std::cout<<this->InterMediate[sample][0][i];
  //		}
  //		std::cout<<" ";
  //	};
  //	std::cout<<"\n";
  for (int i = 0; i < alleleType1.size(); i++) {
    for (int j = 0; j < alleleType2.size(); j++) {
      curHap.reset(new short[2]);
      curHap[0] = alleleType1[i];
      curHap[1] = alleleType2[j];
      OneSampleHaplotypes.push_back(curHap);
    }
  }
  std::vector<std::vector<int> > OneSampleCombination;
  OneSampleGetCombination(OneSampleHaplotypes, this->data->getNumOfChrSet(),
                          OneSampleCombination);

  //	std::vector< boost::shared_ptr<short[]> > v;
  //	std::vector<int> repeatHap;
  //	for(int i=0;i<this->data->getNumOfChrSet();i++){
  //		for(int j=0;j<this->data->getNumOfChrSet();j++){
  //			boost::shared_ptr<short[]> hp(new short[2]);
  //			hp[0]=this->InterMediate[sample][0][i];
  //			hp[1]=this->InterMediate[sample][1][j];
  //			std::cout<<"hap:"<<hp[0]<<hp[1]<<"\n";
  //			if(ExistsHap(hp,v)){
  //				repeatHap.push_back(v.size());
  //			}
  //			v.push_back(hp);
  //		}
  //	}
  //	std::cout<<"all v:";
  //	for(int i=0;i<v.size();i++){
  //		std::cout<<v[i][0]<<v[i][1]<<",";
  //	}
  //	std::cout<<"\n";
  //	std::cout<<"combination.size:"<<this->combination.size()<<"\n";
  for (int i = 0; i < OneSampleCombination.size(); i++) {
    std::vector<int> cur = OneSampleCombination[i];
    std::vector<short> snp1;
    std::vector<short> snp2;
    for (int k = 0; k < cur.size(); k++) {
      int curIdx = cur[k];
      snp1.push_back(OneSampleHaplotypes[curIdx][0]);
      snp2.push_back(OneSampleHaplotypes[curIdx][1]);
    }
    std::sort(snp1.begin(), snp1.end());
    std::sort(snp2.begin(), snp2.end());
    //		std::cout<<"snp1:";
    //		for(int s=0;s<snp1.size();s++){
    //			std::cout<<snp1[s]<<"/";
    //		};
    //		std::cout<<"\nsnp2:";
    //		for(int s=0;s<snp2.size();s++){
    //			std::cout<<snp2[s]<<"/";
    //		}
    //		std::cout<<"\n";
    //		std::cout<<"Genotype:";
    //		for(int p=0;p<this->data->getNumOfChrSet();p++){
    //			std::cout<<this->InterMediate[sample][0][p]<<"/";
    //		}
    //		std::cout<<" ";
    //		for(int p=0;p<this->data->getNumOfChrSet();p++){
    //			std::cout<<this->InterMediate[sample][1][p]<<"/";
    //		}
    //		std::cout<<"\n";

    //		std::cout<<"all haplotypes:";
    //		for(int h=0;h<OneGenotypeExpandedHaplo::haploType.size();h++){
    //			std::cout<<OneGenotypeExpandedHaplo::haploType[h][0]<<OneGenotypeExpandedHaplo::haploType[h][1]<<",";
    //		}
    //		std::cout<<"\n";

    if (GenotypeEqual(snp1, this->InterMediate, sample, 0,
                      this->data->getNumOfChrSet()) &&
        GenotypeEqual(snp2, this->InterMediate, sample, 1,
                      this->data->getNumOfChrSet())) {
      HaploCombination hc;
      for (int n = 0; n < cur.size(); n++) {
        int curIdx = cur[n];
        int hap = getHaploIdx(OneSampleHaplotypes[curIdx],
                              OneGenotypeExpandedHaplo::haploType);
        if (hap == -1) {
          int a = 0;
        }
        hc.hapIdx.push_back(hap);
      }
      if (!ExistsHaploCombination(hc, res.hp)) {
        res.hp.push_back(hc);

        //				std::cout<<"hapIdx combination:";
        //				for(int
        //whichhc=0;whichhc<hc.hapIdx.size();whichhc++){
        //						std::cout<<hc.hapIdx[whichhc]<<",";
        //				};
        //				std::cout<<"\n";

        //				std::cout<<"***Add Combination:";
        //				for(int
        //whichhc=0;whichhc<hc.hapIdx.size();whichhc++){
        //					std::cout<<OneGenotypeExpandedHaplo::haploType[hc.hapIdx[whichhc]][0]<<OneGenotypeExpandedHaplo::haploType[hc.hapIdx[whichhc]][1]<<",";
        //				}
        //				std::cout<<"\n";
      };
    }
    // else
    //		{
    //			std::cout<<"&&&Drop Combination:\n";
    //		}
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

// OneGenotypeExpandedHaplo HaplotypeEM::OneGenoExpandHaplo(int sample) {
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
  for (int i = 0; i < s.hapIdx.size(); i++) {
    if (s.hapIdx[i] == hapidx) return true;
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
    std::vector<double> M(this->UniqueGenotypeCount.size(), 0);
    //  boost::shared_ptr<double[]> M(new
    // double[this->UniqueGenotypeCount.size()]);
    for (int i = 0; i < this->UniqueGenotypeCount.size(); i++) {
      for (int j = 0; j < this->Expanded[i].hp.size(); j++) {
        double m = 1;
        for (int k = 0; k < this->data->getNumOfChrSet(); k++) {
          m *= H[this->Expanded[i].hp[j].hapIdx[k]];
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
            double u = 1;
            for (int kk = 0; kk < this->data->getNumOfChrSet(); kk++) {
              u *= H[this->Expanded[j].hp[k].hapIdx[kk]];
            }
            U += u;
          }
        }
        G += (double)this->UniqueGenotypeCount[j] /
             (double)this->data->getSampleNum() /
             (double)this->data->getNumOfChrSet() * (U / M[j]);
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
  //  for(int i=0;i<OneGenotypeExpandedHaplo::hapfreq.size();i++){
  //	 std::cout<<"Haplotype:"<< OneGenotypeExpandedHaplo::haploType[i][0]<<
  //OneGenotypeExpandedHaplo::haploType[i][1]<<",freq:"
  //			 <<OneGenotypeExpandedHaplo::hapfreq[i]<<"\n";
  //  }
}

void HaplotypeEM::getFinalHap() {
  for (int i = 0; i < this->UniqueGenotypeCount.size(); i++) {
    this->Expanded[i].freq = 0;
    for (int j = 0; j < this->Expanded[i].hp.size(); j++) {
      double freq = 1;
      for (int k = 0; k < this->data->getNumOfChrSet(); k++) {
        int h = this->Expanded[i].hp[j].hapIdx[k];
        freq *= OneGenotypeExpandedHaplo::hapfreq[h];
      };
      if (this->Expanded[i].freq < freq) {
        this->Expanded[i].finalHapIdxInHp = j;
        this->Expanded[i].freq = freq;
      }
    }
  }
}

int locate(boost::multi_array<short, 3> inter, int sample, int ploidy,
           short allele, std::vector<bool>& isPhased) {
  for (int p = 0; p < ploidy; p++) {
    if (inter[sample][0][p] == allele && !isPhased[p]) {
      isPhased[p] = true;
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
    for (int sampleIdx = 0; sampleIdx < samples.size(); sampleIdx++) {
      int finalHapIdx = this->Expanded[i].finalHapIdxInHp;
      HaploCombination CurCombination = this->Expanded[i].hp[finalHapIdx];
      std::vector<bool> isPhased;
      isPhased.resize(this->data->getNumOfChrSet(), false);
      //    	for(int pp=0;pp<2;pp++){
      //			for(int p=0;p<this->data->getNumOfChrSet();p++){
      //				std::cout<<this->InterMediate[samples[sampleIdx]][pp][p];
      //			}
      //			std::cout<<" ";
      //    	};
      //    	std::cout<<"\n";
      for (int p = 0; p < this->data->getNumOfChrSet(); p++) {
        short CurHapAllele1 =
            OneGenotypeExpandedHaplo::haploType[CurCombination.hapIdx[p]][0];
        short CurHapAllele2 =
            OneGenotypeExpandedHaplo::haploType[CurCombination.hapIdx[p]][1];
        int whichPloidy =
            locate(this->InterMediate, samples[sampleIdx],
                   this->data->getNumOfChrSet(), CurHapAllele1, isPhased);
        //    		std::cout<<"a1:"<<CurHapAllele1<<",a2:"<<CurHapAllele2<<",p:"<<whichPloidy<<"\n";
        BOOST_ASSERT(whichPloidy != -1);
        this->PhasedData[samples[sampleIdx]][this->phased][whichPloidy] =
            CurHapAllele2;
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
    //        BOOST_ASSERT((p01 == q01 && p11 == q11) || (p01 == q11 && p11 ==
    // q01));
    //        this->PhasedData[idx][this->phased][0] = p01;
    //        this->PhasedData[idx][this->phased][1] = p11;
    //      } else if (p00 == q10 && p10 == q00) {
    //        BOOST_ASSERT((p01 == q01 && p11 == q11) || (p01 == q11 && p11 ==
    // q01));
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
  std::stringstream p;
  //  std::cout<<"Phased data:\n";
  //  for(int i=0;i<this->data->getSampleNum();i++){
  //	  for(int j=0;j<this->data->getSnpNum();j++){
  //		  for(int k=0;k<this->data->getNumOfChrSet();k++){
  //			  std::cout<<this->PhasedData[i][j][k]<<"/";
  //		  }
  //		  std::cout<<" ";
  //	  }
  //	  std::cout<<"\n";
  //  }
  for (int sample = 0; sample < this->data->getSampleNum(); sample++) {
    if (this->missing[sample]) {
      for (int p = 0; p < this->data->getNumOfChrSet(); p++) {
        this->Results.genotypes[sample][p] = -1;
      }
      continue;
    }

    for (int ploidy = 0; ploidy < this->data->getNumOfChrSet(); ploidy++) {
      p.str("");
      for (int snp = 0; snp < this->SnpIdx.size(); snp++) {
        p << this->PhasedData[sample][snp][ploidy] << ",";
      }
      if (hm.find(p.str()) == hm.end()) {
        haplo.reset(new short[this->SnpIdx.size()]);
        for (int snp = 0; snp < this->SnpIdx.size(); snp++) {
          haplo[snp] = this->PhasedData[sample][snp][ploidy];
        }
        this->Results.haplotypes.push_back(haplo);
        hm[p.str()] = idx++;
      }
      this->Results.genotypes[sample][ploidy] = hm[p.str()];
    }
  };
  //  for(int i=0;i<this->data->getSampleNum();i++){
  //	  	 std::cout<<"sample"<<i;
  //	  for(int j=0;j<this->data->getNumOfChrSet();j++){
  //		  std::cout<<" "<<this->Results.genotypes[i][j];
  //	  }
  //	  std::cout<<"\n";
  //  }

  int haploNum = hm.size();
  this->Results.CaseCount.reset(new int[haploNum]);
  this->Results.ControlCount.reset(new int[haploNum]);
  this->Results.BothCount.reset(new int[haploNum]);
  for (int i = 0; i < this->Results.haplotypes.size(); i++) {
    this->Results.CaseCount[i] = 0;
    this->Results.ControlCount[i] = 0;
    this->Results.BothCount[i] = 0;
  }

  for (int sample = 0; sample < this->data->getSampleNum(); sample++) {
    if (this->missing[sample]) continue;
    if (this->data->vQuantitativeTrait.size() == 0) {
      if (CASE == this->data->vLabel[sample]) {
        for (int ploidy = 0; ploidy < this->data->getNumOfChrSet(); ploidy++) {
          int h = this->Results.genotypes[sample][ploidy];
          this->Results.CaseCount[h]++;
          this->Results.BothCount[h]++;
        }
      } else if (CONTROL == this->data->vLabel[sample]) {
        for (int ploidy = 0; ploidy < this->data->getNumOfChrSet(); ploidy++) {
          int h = this->Results.genotypes[sample][ploidy];
          this->Results.ControlCount[h]++;
          this->Results.BothCount[h]++;
        }
      }
    } else {
      for (int ploidy = 0; ploidy < this->data->getNumOfChrSet(); ploidy++) {
        int h = this->Results.genotypes[sample][ploidy];
        this->Results.BothCount[h]++;
      }
    }
  };

  //	std::cout<<"\nHaplotypes:\n";
  //	for(int k=0;k<this->Results.haplotypes.size();k++){
  //		for(int i=0;i<this->SnpIdx.size();i++){
  //	    	std::cout<<this->Results.haplotypes[k][i];
  //		}
  //		std::cout<<":"<< (this->Results.BothCount[k])<<"\n";
  //	}
}

void HaplotypeEM::startHaplotypeAnalysis() {
  // this->sortGenotype();
  //	std::cout<<"generating combination...\n";
  //	this->getCombination();
  if (0 != this->NonmissingCase + this->NonmissingCtrl) {
    while (this->phased < this->SnpIdx.size()) {
      if (!this->silent && this->phased % 10 == 0) {
        int per = 100 * (double)this->phased / (double)this->SnpIdx.size();
        printf("\rProgress:%d%%", per);
        fflush(stdout);
      }
      this->GenerateInterMediate();
      // this->sortInterMediate();
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
  }
  this->getResults();
}

} /* namespace SHEsis */
