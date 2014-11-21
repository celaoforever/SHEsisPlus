/*
 * HaplotypeDiploid.cpp
 *
 *  Created on: Sep 9, 2014
 *      Author: ada
 */

#include "HaplotypeDiploid.h"
#include "fisher.h"
#include "utility.h"
namespace SHEsis {
std::vector<boost::shared_ptr<short[]> > OneGenotypeExpandedHaplo::haploType;
std::vector<double> OneGenotypeExpandedHaplo::hapfreq;
HaplotypeDiploid::HaplotypeDiploid(boost::shared_ptr<SHEsisData> data)
    : HaplotypeBase(data),
      phased(1),
      PhasedData(boost::extents[data->getSampleNum()][data->getSnpNum()][2]),
      InterMediate(boost::extents[data->getSampleNum()][2][2]),
      err(0.0001),
      missing(new bool[data->getSampleNum()]) {
  BOOST_ASSERT(data->getNumOfChrSet() == 2);
  for (int i = 0; i < this->data->getSnpNum(); i++) {
    this->SnpIdx.push_back(i);
  }
  for (int i = 0; i < data->getSampleNum(); i++) {
    for (int j = 0; j < 2; j++) {
      PhasedData[i][0][j] = data->mGenotype[i][this->SnpIdx[0]][j];
    };
  };
  this->statMissing();
}

HaplotypeDiploid::HaplotypeDiploid(boost::shared_ptr<SHEsisData> data, int Snp,
                                   std::vector<short> mask)
    : HaplotypeBase(data, mask),
      phased(1),
      PhasedData(boost::extents[data->getSampleNum()][data->getSnpNum()][2]),
      InterMediate(boost::extents[data->getSampleNum()][2][2]),
      err(0.00001),
      missing(new bool[data->getSampleNum()]) {
  BOOST_ASSERT(data->getNumOfChrSet() == 2);

  for (int i = 0; i < this->mask.size(); i++) {
    if (this->mask[i]) {
      this->SnpIdx.push_back(i);
    };
  };
  BOOST_ASSERT(this->SnpIdx.size() == Snp);

  for (int i = 0; i < data->getSampleNum(); i++) {
    for (int j = 0; j < 2; j++) {
      PhasedData[i][0][j] = data->mGenotype[i][this->SnpIdx[0]][j];
    };
  };
  this->statMissing();
}

HaplotypeDiploid::~HaplotypeDiploid() {}

void HaplotypeDiploid::statMissing() {
  for (int sample = 0; sample < this->data->getSampleNum(); sample++) {
    this->missing[sample] = false;
    for (int snp = 0; snp < this->SnpIdx.size(); snp++) {
      if (0 == this->data->mGenotype[sample][this->SnpIdx[snp]][0] ||
          0 == this->data->mGenotype[sample][this->SnpIdx[snp]][1])
        this->missing[sample] = true;
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
                   int snp) {
  if (((data[sample1][snp][0] == data[sample2][snp][0]) &&
       (data[sample1][snp][1] == data[sample2][snp][1])) ||
      ((data[sample1][snp][1] == data[sample2][snp][0]) &&
       (data[sample1][snp][0] == data[sample2][snp][1]))) {
    return true;
  } else
    return false;
}

int getPhenotypeCode(std::vector<std::string>& v, std::string str) {
  for (int i = 0; i < v.size(); i++) {
    if (strcmp(v[i].c_str(), str.c_str()) == 0) return i;
  }
  v.push_back(str);
  return (v.size() - 1);
}
void HaplotypeDiploid::ReturnGenotypeCode(int sample, short& geno1,
                                          short& geno2) {
  std::stringstream p1("");
  std::stringstream p2("");
  for (int i = 0; i < this->phased; i++) {
    p1 << this->PhasedData[sample][i][0];
  };
  geno1 = getPhenotypeCode(this->InterMediateGenoCode, p1.str());

  for (int i = 0; i < this->phased; i++) {
    p2 << this->PhasedData[sample][i][1];
  };
  geno2 = getPhenotypeCode(this->InterMediateGenoCode, p2.str());
}

void HaplotypeDiploid::GenerateInterMediate() {
  this->InterMediateGenoCode.clear();
  InterMediateGenoCode.push_back("place holder");
  for (int i = 0; i < data->getSampleNum(); i++) {
    short g1, g2;
    if (missing[i]) {
      g1 = g2 = 0;
    } else {
      this->ReturnGenotypeCode(i, g1, g2);
    }
    this->InterMediate[i][0][0] = g1;
    this->InterMediate[i][0][1] = g2;
    this->InterMediate[i][1][0] =
        this->data->mGenotype[i][this->SnpIdx[this->phased]][0];
    this->InterMediate[i][1][1] =
        this->data->mGenotype[i][this->SnpIdx[this->phased]][1];
  }
}

void HaplotypeDiploid::GenerateUniqueGenotype() {
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
                        0) &&
          GenotypeEqual(this->InterMediate, idx, this->UniqueGenotypeIdx[i],
                        1)) {
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
                                 std::vector<int> sampleIdx, int snp) {
  std::vector<short> alleleType;
  boost::unordered_map<short, int> hm;
  for (int k = 0; k < sampleIdx.size(); k++) {
    for (int i = 0; i < 2; i++) {
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

bool ExistsHaploPair(HaploPair& ap, std::vector<HaploPair>& expanded) {
  for (int i = 0; i < expanded.size(); i++) {
    if (ap == expanded[i]) return true;
  }
  return false;
}

void HaplotypeDiploid::generateAllPossibleHap() {
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

OneGenotypeExpandedHaplo HaplotypeDiploid::OneGenoExpandHaplo(int sample) {
  OneGenotypeExpandedHaplo OneGenoHp;
  for (int i = 0; i < OneGenoHp.haploType.size(); i++) {
    if (compatitable(this->InterMediate, sample,
                     OneGenotypeExpandedHaplo::haploType[i])) {
      HaploPair hp;
      hp.hap1 = i;
      hp.hap2 = PickTheOtherHaplo(this->InterMediate, sample,
                                  OneGenotypeExpandedHaplo::haploType, i);
      BOOST_ASSERT(-1 != hp.hap2);
      if (!ExistsHaploPair(hp, OneGenoHp.hp)) OneGenoHp.hp.push_back(hp);
    };
  }
  return OneGenoHp;
}

void HaplotypeDiploid::ExpandAllGenotype() {
  this->Expanded.clear();
  for (int i = 0; i < this->UniqueGenotypeIdx.size(); i++) {
    this->Expanded.push_back(OneGenoExpandHaplo(this->UniqueGenotypeIdx[i]));
  }
}
bool containsHap(HaploPair s, int hapidx) {
  return (s.hap1 == hapidx || s.hap2 == hapidx);
}

void HaplotypeDiploid::CalculateFreq() {
  int hapcount = OneGenotypeExpandedHaplo::haploType.size();
  OneGenotypeExpandedHaplo::hapfreq.clear();
  // std::cout<<"hapfreq.size():"<<OneGenotypeExpandedHaplo::hapfreq.size()<<"\n";
  for (int i = 0; i < hapcount; i++) {
    OneGenotypeExpandedHaplo::hapfreq.push_back(1.0 / (double)hapcount);
  };

  std::vector<double> H(OneGenotypeExpandedHaplo::haploType.size(),
                        1.0 / (double)hapcount);

  double E = 1;
  while (E > this->err) {
    std::vector<double> M(this->UniqueGenotypeCount.size(), 0);
    for (int i = 0; i < this->UniqueGenotypeCount.size(); i++) {
      for (int j = 0; j < this->Expanded[i].hp.size(); j++) {
        int h1 = this->Expanded[i].hp[j].hap1;
        int h2 = this->Expanded[i].hp[j].hap2;
        M[i] += H[h1] * H[h2];
      }
    }
    for (int n = 0; n < hapcount; n++) {
      double G = 0;
      for (int j = 0; j < this->UniqueGenotypeCount.size(); j++) {
        double U = 0;
        for (int k = 0; k < this->Expanded[j].hp.size(); k++) {
          if (containsHap(this->Expanded[j].hp[k], n)) {
            U += H[this->Expanded[j].hp[k].hap1] *
                 H[this->Expanded[j].hp[k].hap2];
          }
        }
        G += (double)this->UniqueGenotypeCount[j] /
             (double)this->data->getSampleNum() / 2.0 * (U / M[j]);
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

void HaplotypeDiploid::getFinalHap() {
  for (int i = 0; i < this->UniqueGenotypeCount.size(); i++) {
    this->Expanded[i].freq = 0;
    for (int j = 0; j < this->Expanded[i].hp.size(); j++) {
      int h1 = this->Expanded[i].hp[j].hap1;
      int h2 = this->Expanded[i].hp[j].hap2;
      double freq = OneGenotypeExpandedHaplo::hapfreq[h1] *
                    OneGenotypeExpandedHaplo::hapfreq[h2];
      if (this->Expanded[i].freq < freq) {
        this->Expanded[i].finalhap1 = h1;
        this->Expanded[i].finalhap2 = h2;
        this->Expanded[i].freq = freq;
      }
    }
  }
}

std::vector<int> HaplotypeDiploid::getSampleIdx(int genotype) {
  std::vector<int> res;
  for (int i = 0; i < this->data->getSampleNum(); i++) {
    if (this->Sample2Genotype[i] == genotype) res.push_back(i);
  }
  return res;
}
void HaplotypeDiploid::PhaseCurrent() {
  for (int i = 0; i < this->UniqueGenotypeCount.size(); i++) {
    std::vector<int> samples = getSampleIdx(i);
    int h1 = this->Expanded[i].finalhap1;
    int h2 = this->Expanded[i].finalhap2;
    short p00 = OneGenotypeExpandedHaplo::haploType[h1][0];
    short p01 = OneGenotypeExpandedHaplo::haploType[h1][1];
    short p10 = OneGenotypeExpandedHaplo::haploType[h2][0];
    short p11 = OneGenotypeExpandedHaplo::haploType[h2][1];
    for (int j = 0; j < samples.size(); j++) {
      int idx = samples[j];
      short q00 = this->InterMediate[idx][0][0];
      short q01 = this->InterMediate[idx][1][0];
      short q10 = this->InterMediate[idx][0][1];
      short q11 = this->InterMediate[idx][1][1];
      if (p00 == q00 && p10 == q10) {
        BOOST_ASSERT((p01 == q01 && p11 == q11) || (p01 == q11 && p11 == q01));
        this->PhasedData[idx][this->phased][0] = p01;
        this->PhasedData[idx][this->phased][1] = p11;
      } else if (p00 == q10 && p10 == q00) {
        BOOST_ASSERT((p01 == q01 && p11 == q11) || (p01 == q11 && p11 == q01));
        this->PhasedData[idx][this->phased][0] = p11;
        this->PhasedData[idx][this->phased][1] = p01;

      } else {
        BOOST_ASSERT(1 == 0);
      }
    }
  }
  this->phased++;
};

void HaplotypeDiploid::getResults() {
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
    } else {
      for (int ploidy = 0; ploidy < this->data->getNumOfChrSet(); ploidy++) {
        int h = this->Results.genotypes[sample][ploidy];
        this->Results.BothCount[h]++;
      }
    }
  };
  std::cout << "\nHaplotypes:\n";
  for (int k = 0; k < this->Results.haplotypes.size(); k++) {
    for (int i = 0; i < this->SnpIdx.size(); i++) {
      std::cout << this->Results.haplotypes[k][i];
    }
    std::cout << ":"
              << (this->Results.ControlCount[k] + this->Results.CaseCount[k])
              << "\n";
  }
}

void HaplotypeDiploid::startHaplotypeAnalysis() {
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
