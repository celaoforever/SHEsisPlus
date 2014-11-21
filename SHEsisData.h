/*
 * SHEsisData.h
 *
 *  Created on: Aug 3, 2014
 *      Author: ionadmin
 */

#ifndef SHESISDATA_H_
#define SHESISDATA_H_
#include <boost/multi_array.hpp>
#include <boost/unordered_map.hpp>
#include <boost/assert.hpp>
#include <vector>
#include <string>

#define GENOTYPE_MISSING 0
namespace SHEsis {

struct Missing {
  Missing() : CaseAlleleNum(0), CtrlAlleleNum(0) {}
  int CaseAlleleNum;
  int CtrlAlleleNum;
};

struct LocusInfo {
  LocusInfo() : AlleleCallrate(0), GenoCallrate(0) {};
  boost::unordered_map<short, double> CaseAlleleCount;
  boost::unordered_map<short, double> ControlAlleleCount;
  boost::unordered_map<std::string, double> CaseGenotypeCount;
  boost::unordered_map<std::string, double> ControlGenotypeCount;

  boost::unordered_map<short, double> BothAlleleCount;
  boost::unordered_map<std::string, double> BothGenotypeCount;
  double AlleleCallrate;
  double GenoCallrate;
  int getAlleleIndex(short a) {
    BOOST_ASSERT(BothAlleleCount.end() != BothAlleleCount.find(a));
    boost::unordered_map<short, double>::iterator iter;
    int idx = 0;
    for (iter = BothAlleleCount.begin(); iter != BothAlleleCount.end();
         iter++) {
      if (a == iter->first) return idx;
      idx++;
    }
    return -1;
  };
  short getAlleleType(int idx) {
    BOOST_ASSERT(idx < BothAlleleCount.size());
    boost::unordered_map<short, double>::iterator iter;
    int i = 0;
    for (iter = BothAlleleCount.begin(); iter != BothAlleleCount.end();
         iter++) {
      if (i == idx) return iter->first;
      i++;
    }
    return -1;
  }
};

typedef enum { MISSING, CONTROL, CASE } SampleStatus;

class SHEsisData {
 public:
  SHEsisData(int SampleNum, int SnpNum, int NumOfChrSet);
  int getNumOfChrSet() {
    return this->NumOfChrSet;
  };
  int getSampleNum() {
    return this->SampleNum;
  };
  int getSnpNum() {
    return this->SnpNum;
  };
  virtual ~SHEsisData();
  boost::multi_array<short, 3> mGenotype;
  std::vector<double> vQuantitativeTrait;
  std::vector<SampleStatus> vLabel;
  std::vector<LocusInfo> vLocusInfo;
  std::vector<std::string> vLocusName;
  void statCount(std::vector<SampleStatus>& label);
  void statCount();  // for qtl
  void printLocusInfo();
  int getCaseNum();
  int getControlNum();
  short GetAlleleCode(std::string const val);
  void setLocusName(int snpidx, std::string s);
  double getCaseCallrate(int snp);
  double getControlCallrate(int snp);
  double getCallrate(int snp);
  std::string getStrFromSortedGenotype(std::vector<short> v);
  std::string getallele(short a) {
    if (code2allele.find(a) != code2allele.end()) return code2allele[a];
    return "NA";
  };
  std::string getOriginGenotype(std::string geno);

 protected:
  boost::unordered_map<short, std::string> code2allele;
  int codeIdx;
  void getCaseAndControlNum();
  int CaseNum;
  int ControlNum;
  std::vector<Missing> missingalleles;
  const int SampleNum;
  const int SnpNum;
  const int NumOfChrSet;
};
};

#endif /* SHESISDATA_H_ */
