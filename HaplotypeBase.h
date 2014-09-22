/*
 * HaplotypeBase.h
 *
 *  Created on: Sep 9, 2014
 *      Author: ada
 */

#ifndef HAPLOTYPEBASE_H_
#define HAPLOTYPEBASE_H_
#include <boost/shared_ptr.hpp>
#include "SHEsisData.h"
#include "fisher.h"
#include "utility.h"
namespace SHEsis {

struct singleHapRes {
  singleHapRes()
      : chisquare(-1), fisherp(-1), pearsonp(-1), OR(-1), orlow(-1), orUp(-1) {}
  double chisquare;
  double fisherp;
  double pearsonp;
  double OR;
  double orlow;
  double orUp;
};
struct HapTestResult {
  std::vector<boost::shared_ptr<short[]> > haplotypes;
  std::vector<singleHapRes> singleHap;
  boost::multi_array<short, 2> genotypes;
  boost::shared_ptr<int[]> CaseCount;
  boost::shared_ptr<int[]> ControlCount;
  double FisherP;
  double PearsonP;
  double ChiSquare;
  HapTestResult(int sample, int ploidy)
      : genotypes(boost::extents[sample][ploidy]),
        FisherP(-1),
        PearsonP(-1),
        ChiSquare(-1) {};
  ~HapTestResult() { haplotypes.clear(); }
};

class HaplotypeBase {
 public:
  boost::shared_ptr<SHEsisData> data;
  HaplotypeBase(boost::shared_ptr<SHEsisData> data, std::vector<short> mask)
      : data(data),
        mask(mask),
        Results(data->getSampleNum(), data->getNumOfChrSet()),
        silent(true),
        freqthreshold(0.03) {};
  HaplotypeBase(boost::shared_ptr<SHEsisData> data)
      : data(data),
        Results(data->getSampleNum(), data->getNumOfChrSet()),
        silent(true),
        freqthreshold(0.03) {};
  virtual ~HaplotypeBase() {
    this->mask.clear();
    this->SnpIdx.clear();
  }
  virtual void startHaplotypeAnalysis() {};
  void AssociationTest();
  void setSilent(bool b) {
    this->silent = b;
  };
  void setFreqThreshold(double t) {
    this->freqthreshold = t;
  };
  std::string reporthtml();
  std::string reporthtmltable();
  HapTestResult Results;

 protected:
  std::vector<short> mask;
  std::vector<int> SnpIdx;
  bool silent;
  double freqthreshold;
};

} /* namespace SHEsis */

#endif /* HAPLOTYPEBASE_H_ */
