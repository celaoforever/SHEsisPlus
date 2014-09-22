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
  double AllelePearsonP;
  double AlleleFisherP;
  double AlleleChiSquare;
  double AlleleOddsRatio;
  double AlleleOddsRatioLowLimit;
  double AlleleOddsRatioUpLimit;
  double AllelePermutationP;
};

class AssociationTest {
 public:
  boost::shared_ptr<SHEsisData> data;
  std::vector<LocusAssiciationTestResult> vAssocationTestResult;
  AssociationTest(boost::shared_ptr<SHEsisData> data);
  void setPermutationTimes(int p) {
    this->NumOfPermutation = p;
  };
  virtual ~AssociationTest();
  void printAssociationTestResults();
  void report();
  void permutation();
  void association();
  std::string reporthtml();

 private:
  int NumOfPermutation;
  std::vector<double> PermutationPAllele;
  std::vector<double> PermutationPGenotype;
  std::vector<SampleStatus> vPermutateLabel;
  void SingleSnpTestAllele(int iSnp, double& FisherP, double& PearsonP,
                           double& ChiSquare, double& oddsRatio,
                           double& ORLowLimit, double& ORUpLimit);
  void SingleSnpTestGenotype(int iSnp, double& FisherP, double& PearsonP,
                             double& ChiSquare);
  void AssociationTestForAllSnpsAllele();
  void AssociationTestForAllSnpsGenotype();
  std::string reporthtmlAllele();
  std::string reporthtmlGenotype();
};

} /* namespace SHEsis */

#endif /* ASSOCIATIONTEST_H_ */
