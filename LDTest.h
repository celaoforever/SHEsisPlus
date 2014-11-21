/*
 * LDTest.h
 *
 *  Created on: Sep 17, 2014
 *      Author: ada
 */

#ifndef LDTEST_H_
#define LDTEST_H_
#include "HaplotypeBase.h"
#include <boost/shared_ptr.hpp>
#include "BMP.h"
namespace SHEsis {

typedef enum { LD_IN_CASE, LD_IN_CTRL, LD_IN_BOTH } LD_TYPE;
struct NonmissingSampleAlleleCount {
  boost::unordered_map<short, int> casecount;
  boost::unordered_map<short, int> ctrlcount;
};

class LDTest {
 public:
  LDTest(boost::shared_ptr<SHEsisData> data, std::string path);
  virtual ~LDTest();
  void AllLociLDtest();
  void DrawLDMapDandR2();
  void printRes();
  void setLDType(LD_TYPE t) {
    this->ldtype = t;
  };
  void setForceSAT(bool b) {
    this->bForceSAT = b;
  };
  std::string reporthtml();
  std::string reporttxt();

 private:
  void TwoLociLDTest(int snp1, int snp2, LD_TYPE type, double& R2, double& D);
  void DrawLDMap(int type);
  void StatAllelesInNonmissingSample(boost::shared_ptr<bool[]> missing,
                                     int snp1, int snp2);
  void resetMissing();
  NonmissingSampleAlleleCount Snp1AlleleCount;
  NonmissingSampleAlleleCount Snp2AlleleCount;
  boost::shared_ptr<SHEsisData> data;
  boost::shared_ptr<HaplotypeBase> hp;
  boost::multi_array<double, 2> resD;
  boost::multi_array<double, 2> resR2;
  BMP* ldmap;
  std::string path;
  LD_TYPE ldtype;
  bool bForceSAT;
};

} /* namespace SHEsis */

#endif /* LDTEST_H_ */
