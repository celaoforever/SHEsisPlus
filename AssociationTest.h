/*
 * AssociationTest.h
 *
 *  Created on: Aug 8, 2014
 *      Author: ada
 */

#ifndef ASSOCIATIONTEST_H_
#define ASSOCIATIONTEST_H_
#include "SHEsisData.h"
#include <boost/shared_ptr.hpp>
namespace SHEsis {

struct LocusAssiciationTestResult {
  double GenoTypePearsonP;
  double GenotypeFisherP;
  double GenotypeChiSquare;
  double GenotypePermutationP;  // pearson's p
  double GenotypeHolmP;
  double GenotypeSidakSSP;
  double GenotypeSidakSDP;
  double GenotypeBHP;
  double GenotypeBYP;
  double AllelePearsonP;
  double AlleleFisherP;
  double AlleleChiSquare;
  double AlleleOddsRatio;
  double AlleleOddsRatioLowLimit;
  double AlleleOddsRatioUpLimit;
  double AllelePermutationP;
  double AlleleHolmP;
  double AlleleSidakSSP;
  double AlleleSidakSDP;
  double AlleleBHP;
  double AlleleBYP;
};
// void HolmCorrection(std::vector<double>& p, std::vector<double>& adjusted);
// void SidakSSCorrection(std::vector<double>& p, std::vector<double>&
// adjusted);
// void SidakSDCorrection(std::vector<double>& p, std::vector<double>&
// adjusted);
// void BHCorrection(std::vector<double>& p, std::vector<double>& adjusted);
// void BYCorrection(std::vector<double>& p, std::vector<double>& adjusted);
class AssociationTest {
 public:
  boost::shared_ptr<SHEsisData> data;
  std::vector<LocusAssiciationTestResult> vAssocationTestResult;
  AssociationTest(boost::shared_ptr<SHEsisData> data);
  void setPermutationTimes(int p) {
    this->NumOfPermutation = p;
  };
  void setAdjust(bool b) {
    this->adjust = b;
  };
  virtual ~AssociationTest();
  void printAssociationTestResults();
  void permutation();
  void association();
  std::string reporthtml();
  std::string reporttxt();

 private:
  int NumOfPermutation;
  std::vector<double> PermutationPAllele;
  std::vector<double> PermutationPGenotype;
  std::vector<SampleStatus> vPermutateLabel;
  bool adjust;
  void SingleSnpTestAllele(int iSnp, double& FisherP, double& PearsonP,
                           double& ChiSquare, double& oddsRatio,
                           double& ORLowLimit, double& ORUpLimit, bool);
  void SingleSnpTestGenotype(int iSnp, double& FisherP, double& PearsonP,
                             double& ChiSquare, bool);
  void AssociationTestForAllSnpsAllele(bool);
  void AssociationTestForAllSnpsGenotype(bool);
  std::string reporthtmlAllele();
  std::string reporthtmlGenotype();
};

} /* namespace SHEsis */

#endif /* ASSOCIATIONTEST_H_ */
