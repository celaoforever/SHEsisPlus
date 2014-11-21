/*
 * QTL.h
 *
 *  Created on: Oct 3, 2014
 *      Author: ada
 */

#ifndef QTL_H_
#define QTL_H_
#include <boost/shared_ptr.hpp>
#include <vector>
#include "SHEsisData.h"
namespace SHEsis {
struct QTLResults {
  QTLResults()
      : p(2),
        ValidSampleNum(0),
        beta(0),
        se(0),
        R2(0),
        T(0),
        allele(0),
        permutatedp(0) {};
  double ValidSampleNum;
  double beta;
  double se;
  double R2;
  double T;
  double p;
  short allele;
  double permutatedp;
  double HolmP;
  double SidakSSP;
  double SidakSDP;
  double BHP;
  double BYP;
  bool operator<(const QTLResults& res) {
    if (p < res.p)
      return true;
    else
      return false;
  };
  bool operator>(const QTLResults& res) {
    if (p > res.p)
      return true;
    else
      return false;
  };
};
class QTL {
 public:
  boost::shared_ptr<SHEsisData> data;
  std::vector<QTLResults> vResults;
  void QTLTest();
  void printRes();
  QTL(boost::shared_ptr<SHEsisData> data);
  virtual ~QTL();
  void setPermutation(int p) {
    this->NumOfPermutation = p;
  };
  void setAdjust(bool b) {
    this->adjust = b;
  };
  void QTLPermutation();
  std::string reporthtml();
  std::string reporttxt();

 private:
  bool adjust;
  int NumOfPermutation;
  std::vector<double> vPermutatedQT;
  QTLResults OneLocusQTLAnalysis(int snp, short allele, std::vector<double> qt);
  int getAlleleCount(int sample, int snp, short allele);
  bool IsGenotypeMissing(int sample, int snp);
  std::vector<short> FindAllele(int snp);
  void QTLTest(std::vector<double> qt);
  std::vector<double> permutationp;
};

} /* namespace SHEsis */

#endif /* QTL_H_ */
