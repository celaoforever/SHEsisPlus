/*
 * HaplotypeDiploid.h
 *
 *  Created on: Sep 9, 2014
 *      Author: ada
 */

#ifndef HAPLOTYPEDEM_H_
#define HAPLOTYPEDEM_H_
#include "HaplotypeBase.h"

namespace SHEsis {

struct HaploCombination {
  std::vector<int> hapIdx;
  bool operator==(const HaploCombination& p) {
    for (int i = 0; i < p.hapIdx.size(); i++) {
      if (p.hapIdx[i] != hapIdx[i]) return false;
    }
    return true;
  }
};

struct OneGenotypeExpandedHaplo {
  static std::vector<boost::shared_ptr<short[]> > haploType;
  static std::vector<double> hapfreq;
  std::vector<HaploCombination> hp;
  double freq;
  int finalHapIdxInHp;
};

class HaplotypeEM : public HaplotypeBase {
 public:
  HaplotypeEM(boost::shared_ptr<SHEsisData> data);
  HaplotypeEM(boost::shared_ptr<SHEsisData> data, int Snp,
              std::vector<short> mask);
  virtual ~HaplotypeEM();
  virtual void startHaplotypeAnalysis();
  void getResults();
  void ShowResults(bool b){this->showResults=b;};
  void setSeed(int s) {this->seed=s;};
  boost::shared_ptr<bool[]> missing;
  boost::multi_array<short, 3> PhasedData;
 private:
  void sortGenotype();
  void sortInterMediate();
  void GenerateUniqueGenotype();
  void GenerateInterMediate();
  void ReturnGenotypeCode(int sample, std::vector<short>& geno);
  OneGenotypeExpandedHaplo OneGenoExpandHaplo(int sample);
  void ExpandAllGenotype();
  void CalculateFreq();
  void generateAllPossibleHap();
  void statMissing();
  void getFinalHap();
  void PhaseCurrent();
  void getCombination();
  OneGenotypeExpandedHaplo oneGenoGetCombination(int sample);
  std::vector<int> getSampleIdx(int genotype);
  std::vector<int> UniqueGenotypeIdx;  // 2 loci, diploid
  std::vector<int> UniqueGenotypeCount;
  boost::multi_array<short, 3> InterMediate;
  std::vector<std::string> InterMediateGenoCode;
  std::vector<int> CurGenotypeCount;
  boost::shared_ptr<int[]> Sample2Genotype;
  std::vector<OneGenotypeExpandedHaplo> Expanded;
  std::vector<std::vector<int> > combination;
  std::vector<double> getInitialFreq(int hapcount,int seed);
  int seed;
  int phased;
  double err;
  bool showResults;
};

} /* namespace SHEsis */

#endif /* HAPLOTYPEDIPLOID_H_ */
