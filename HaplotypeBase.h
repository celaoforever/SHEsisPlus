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
      : chisquare(-999),
        fisherp(-999),
        pearsonp(-999),
        OR(-999),
        orlow(-999),
        orUp(-999) {}
  double chisquare;
  double fisherp;
  double pearsonp;
  double OR;
  double orlow;
  double orUp;
  double HolmP;
  double SidakSSP;
  double SidakSDP;
  double BHP;
  double BYP;
};
struct singHapQTLRes {
  singHapQTLRes() : p(1), ValidSampleNum(0), beta(0), se(0), R2(0), T(0) {};
  double ValidSampleNum;
  double beta;
  double se;
  double R2;
  double T;
  double p;
  double HolmP;
  double SidakSSP;
  double SidakSDP;
  double BHP;
  double BYP;
};
struct HapTestResult {
  std::vector<boost::shared_ptr<short[]> > haplotypes;
  std::vector<singleHapRes> singleHap;
  std::vector<singHapQTLRes> singleHapQTL;
  boost::multi_array<int, 2> genotypes;
  boost::shared_ptr<int[]> CaseCount;
  boost::shared_ptr<int[]> ControlCount;
  boost::shared_ptr<int[]> BothCount;
  double FisherP;
  double PearsonP;
  double ChiSquare;
  HapTestResult(int sample, int ploidy)
      : genotypes(boost::extents[sample][ploidy]),
        FisherP(-999),
        PearsonP(-999),
        ChiSquare(-999) {};
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
        freqthreshold(0.03),
        adjust(false),
        NonmissingCase(0),
        NonmissingCtrl(0) {};
  HaplotypeBase(boost::shared_ptr<SHEsisData> data)
      : data(data),
        Results(data->getSampleNum(), data->getNumOfChrSet()),
        silent(true),
        freqthreshold(0.03),
        adjust(false),
        NonmissingCase(0),
        NonmissingCtrl(0) {};
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
  void setAdjust(bool b) { this->adjust = b; }

  std::string reporthtml();
  std::string reporttxt();
  HapTestResult Results;
  int NonmissingCase;
  int NonmissingCtrl;

 protected:
  void AssociationTestBinary();
  void AssociationTestQTL();
  singHapQTLRes SingleHaploAssociationTestQTL(int hapIdx);
  int getHaploCount(int sample, int hapIdx);
  bool IsHaploMissing(int sample);

  std::vector<short> mask;
  std::vector<int> SnpIdx;
  std::string reporthtmltableBinary();
  std::string reporthtmltableQTL();
  std::string reporttxttableBinary();
  std::string reporttxttableQTL();
  bool silent;
  bool adjust;
  double freqthreshold;
};

} /* namespace SHEsis */

#endif /* HAPLOTYPEBASE_H_ */
