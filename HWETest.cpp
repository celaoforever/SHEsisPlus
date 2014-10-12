/*
 * HWETest.cpp
 *
 *  Created on: Aug 10, 2014
 *      Author: ionadmin
 */
#include <boost/assert.hpp>
#include "HWETest.h"
#include "Multinominal.h"
#include "fisher.h"
#include "utility.h"
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include "CreatHtmlTable.h"
template std::string convert2string<double>(double v);
namespace SHEsis {

HWETest::HWETest(boost::shared_ptr<SHEsisData> data)
    : data(data), vHWETestResult(data->getSnpNum()) {}

HWETest::~HWETest() {
  this->vHWETestResult.clear();
  this->vCoefficient.clear();
}

void ResetMapValue(boost::unordered_map<short, size_t>& map) {
  boost::unordered_map<short, size_t>::iterator map_it;
  for (map_it = map.begin(); map_it != map.end(); map_it++) {
    map[map_it->first] = 0;
  };
}

std::string getStrFromSortedShort(std::vector<size_t> v) {
  std::stringstream ss;
  for (int i = 0; i < v.size() - 1; i++) {
    ss << v[i] << "/";
  }
  ss << v[v.size() - 1];
  return ss.str();
}

double GetExpectedGenotypeFreq(
    std::string genotype, boost::unordered_map<short, double>& AlleleCount,
    boost::unordered_map<short, size_t>& AlleleType,
    boost::unordered_map<std::string, size_t>& coefficient, int SampleNum,
    int NumOfChrSet) {

  ResetMapValue(AlleleType);
  std::vector<std::string> splitted;
  boost::split(splitted, genotype, boost::is_any_of("/"));
  BOOST_ASSERT(splitted.size() > 0);

  boost::unordered_map<short, size_t>::iterator iter;
  for (int i = 0; i < splitted.size(); i++) {
    iter = AlleleType.find(std::atoi(splitted[i].c_str()));
    BOOST_ASSERT(AlleleType.end() != iter);
    AlleleType[iter->first]++;
  }

  std::vector<size_t> exponent;
  for (iter = AlleleType.begin(); iter != AlleleType.end(); iter++) {
    exponent.push_back(iter->second);
  };
  std::sort(exponent.begin(), exponent.end());
  std::string AlleleTypeStatus = getStrFromSortedShort(exponent);

  boost::unordered_map<std::string, size_t>::iterator iter2 =
      coefficient.find(AlleleTypeStatus);
  size_t CurCoefficient;
  if (coefficient.end() != iter2) {
    CurCoefficient = iter2->second;
  } else {
    CurCoefficient = multi<size_t>(exponent);
    coefficient[AlleleTypeStatus] = CurCoefficient;
  };

  double res = CurCoefficient;
  for (iter = AlleleType.begin(); iter != AlleleType.end(); iter++) {
    short CurAllele = iter->first;
    size_t CurCount = iter->second;
    if (0 == CurCount) continue;
    res *= pow(AlleleCount[CurAllele] / ((double)SampleNum * NumOfChrSet),
               CurCount);
  }

  splitted.clear();
  exponent.clear();

  return res;
}

void HWETest::SingleSnpHWETest(int iSnp, double& CaseChi, double& CasePearsonP,
                               double& CaseFisherP, double& ControlChi,
                               double& ControlPearsonP, double& ControlFisherP,
                               double& BothChi, double& BothPearsonP,
                               double& BothFisherP) {
  BOOST_ASSERT(this->data->vLocusInfo[iSnp].CaseAlleleCount.size() ==
               this->data->vLocusInfo[iSnp].ControlAlleleCount.size());

  boost::unordered_map<short, size_t> AlleleType;
  boost::unordered_map<short, double>::iterator map_it;
  for (map_it = this->data->vLocusInfo[iSnp].CaseAlleleCount.begin();
       map_it != this->data->vLocusInfo[iSnp].CaseAlleleCount.end(); map_it++) {
    AlleleType[map_it->first] = 0;
  };

  // calculate HWE for cases
  boost::unordered_map<std::string, double>::iterator genotype_iter;
  int NumOfRow =
      2;  // 1st row is for observed freq, 2nd row is for expected frequency
  int NumOfCol = this->data->vLocusInfo[iSnp].CaseGenotypeCount.size();
  int totalGenotype = 0;
  double* contigency = new double[NumOfRow * NumOfCol];
  int idx = 0;
  for (genotype_iter = this->data->vLocusInfo[iSnp].CaseGenotypeCount.begin();
       genotype_iter != this->data->vLocusInfo[iSnp].CaseGenotypeCount.end();
       genotype_iter++) {
    double expectedFreq = GetExpectedGenotypeFreq(
        genotype_iter->first, this->data->vLocusInfo[iSnp].CaseAlleleCount,
        AlleleType, this->vCoefficient, this->data->getCaseNum(),
        this->data->getNumOfChrSet());
    BOOST_ASSERT(idx < NumOfRow * NumOfCol);
    contigency[idx++] =
        genotype_iter
            ->second;  //(this->data->getCaseNum()*this->data->getNumOfChrSet());
    BOOST_ASSERT(idx < NumOfRow * NumOfCol);
    contigency[idx++] = expectedFreq;
    totalGenotype += genotype_iter->second;
  };

  // Fisher's exact test:
  CaseChi = 0;
  for (int i = 0; i < NumOfCol * NumOfRow; i = i + 2) {
    if (contigency[i] < 1) contigency[i] *= totalGenotype;
    if (contigency[i + 1] < 1) contigency[i + 1] *= totalGenotype;
    if (contigency[i + 1] != 0)
      CaseChi += (contigency[i + 1] - contigency[i]) *
                 (contigency[i + 1] - contigency[i]) / contigency[i + 1];
  }
  boost::math::chi_squared dist(NumOfCol);
  CasePearsonP = boost::math::cdf(boost::math::complement(dist, CaseChi));

  double expect = -1.0;
  double percnt = 100.0;
  double emin = 0;
  double pre = 0, prt = 0;
  int ws = 300000;
  try {
    fexact(&NumOfRow, &NumOfCol, contigency, &NumOfRow, &expect, &percnt, &emin,
           &prt, &pre, &ws);
    CaseFisherP = pre;
  }
  catch (std::runtime_error&) {
    CaseFisherP = -999;
  }
  // Pearson's ChiSquare test
  // PearsonChiSquareTest(contigency,NumOfRow,NumOfCol,CaseChi,CasePearsonP);
  delete[] contigency;
  contigency = 0;

  // calculate HWE for control
  totalGenotype = 0;
  NumOfRow =
      2;  // 1st row is for observed freq, 2nd row is for expected frequency
  NumOfCol = this->data->vLocusInfo[iSnp].ControlGenotypeCount.size();
  contigency = new double[NumOfRow * NumOfCol];
  idx = 0;
  for (genotype_iter =
           this->data->vLocusInfo[iSnp].ControlGenotypeCount.begin();
       genotype_iter != this->data->vLocusInfo[iSnp].ControlGenotypeCount.end();
       genotype_iter++) {
    double expectedFreq = GetExpectedGenotypeFreq(
        genotype_iter->first, this->data->vLocusInfo[iSnp].ControlAlleleCount,
        AlleleType, this->vCoefficient, this->data->getControlNum(),
        this->data->getNumOfChrSet());
    BOOST_ASSERT(idx < NumOfRow * NumOfCol);
    contigency[idx++] =
        genotype_iter
            ->second;  //(this->data->getControlNum()*this->data->getNumOfChrSet());
    BOOST_ASSERT(idx < NumOfRow * NumOfCol);
    contigency[idx++] = expectedFreq;
    totalGenotype += genotype_iter->second;
  };

  // Fisher's exact test:
  ControlChi = 0;
  for (int i = 0; i < NumOfCol * NumOfRow; i = i + 2) {
    if (contigency[i] < 1) contigency[i] *= totalGenotype;
    if (contigency[i + 1] < 1) contigency[i + 1] *= totalGenotype;
    if (contigency[i + 1] != 0)
      ControlChi += (contigency[i + 1] - contigency[i]) *
                    (contigency[i + 1] - contigency[i]) / contigency[i + 1];
  }
  ControlPearsonP = boost::math::cdf(boost::math::complement(dist, ControlChi));
  try {
    fexact(&NumOfRow, &NumOfCol, contigency, &NumOfRow, &expect, &percnt, &emin,
           &prt, &pre, &ws);
    ControlFisherP = pre;
  }
  catch (std::runtime_error&) {
    ControlFisherP = -999;
  }

  // Pearson's ChiSquare test

  delete[] contigency;
  contigency = 0;

  // calculate HWE for case and control
  totalGenotype = 0;
  NumOfRow =
      2;  // 1st row is for observed freq, 2nd row is for expected frequency
  NumOfCol = this->data->vLocusInfo[iSnp].BothGenotypeCount.size();
  contigency = new double[NumOfRow * NumOfCol];
  idx = 0;
  for (genotype_iter = this->data->vLocusInfo[iSnp].BothGenotypeCount.begin();
       genotype_iter != this->data->vLocusInfo[iSnp].BothGenotypeCount.end();
       genotype_iter++) {
    double expectedFreq = GetExpectedGenotypeFreq(
        genotype_iter->first, this->data->vLocusInfo[iSnp].BothAlleleCount,
        AlleleType, this->vCoefficient, this->data->getSampleNum(),
        this->data->getNumOfChrSet());
    BOOST_ASSERT(idx < NumOfRow * NumOfCol);
    contigency[idx++] = genotype_iter->second;
    BOOST_ASSERT(idx < NumOfRow * NumOfCol);
    contigency[idx++] = expectedFreq;
    totalGenotype += genotype_iter->second;
  };

  // Fisher's exact test:
  BothChi = 0;
  for (int i = 0; i < NumOfCol * NumOfRow; i = i + 2) {
    if (contigency[i] < 1) contigency[i] *= totalGenotype;
    if (contigency[i + 1] < 1) contigency[i + 1] *= totalGenotype;
    if (contigency[i + 1] != 0)
      BothChi += (contigency[i + 1] - contigency[i]) *
                 (contigency[i + 1] - contigency[i]) / contigency[i + 1];
  }
  BothPearsonP = boost::math::cdf(boost::math::complement(dist, BothChi));
  try {
    fexact(&NumOfRow, &NumOfCol, contigency, &NumOfRow, &expect, &percnt, &emin,
           &prt, &pre, &ws);
    BothFisherP = pre;
  }
  catch (std::runtime_error&) {
    BothFisherP = -999;
  }
  // Pearson's ChiSquare test
  delete[] contigency;
  contigency = 0;
}

void HWETest::AllSnpHWETestBinary() {
  if(this->data->vLocusInfo[0].BothAlleleCount.size()==0)
	  this->data->statCount(this->data->vLabel);
  for (int i = 0; i < this->data->getSnpNum(); i++) {
    this->SingleSnpHWETest(i, this->vHWETestResult[i].CaseChiSquare,
                           this->vHWETestResult[i].CasePearsonP,
                           this->vHWETestResult[i].CaseFisherP,
                           this->vHWETestResult[i].ControlChiSquare,
                           this->vHWETestResult[i].ControlPearsonP,
                           this->vHWETestResult[i].ControlFisherP,
                           this->vHWETestResult[i].BothChiSquare,
                           this->vHWETestResult[i].BothPearsonP,
                           this->vHWETestResult[i].BothFisherP);
  };
};

void HWETest::AllSnpHWETestQTL() {
	if(this->data->vLocusInfo[0].BothAlleleCount.size()==0)
		this->data->statCount();
//  this->data->printLocusInfo();
  for (int i = 0; i < this->data->getSnpNum(); i++) {
    this->SingleSnpHWETestBoth(i,this->vHWETestResult[i].BothChiSquare,
                           this->vHWETestResult[i].BothPearsonP,
                           this->vHWETestResult[i].BothFisherP);
  };
};

void HWETest::AllSnpHWETest(){
	if(this->data->vQuantitativeTrait.size()!=0)
		this->AllSnpHWETestQTL();
	else
		this->AllSnpHWETestBinary();
}

void HWETest::SingleSnpHWETestBoth(int iSnp, double& BothChi, double& BothPearsonP,
                               double& BothFisherP) {

  boost::unordered_map<short, size_t> AlleleType;
  boost::unordered_map<short, double>::iterator map_it;
  for (map_it = this->data->vLocusInfo[iSnp].BothAlleleCount.begin();
       map_it != this->data->vLocusInfo[iSnp].BothAlleleCount.end(); map_it++) {
    AlleleType[map_it->first] = 0;
  };

  boost::unordered_map<std::string, double>::iterator genotype_iter;
  int NumOfRow =
      2;  // 1st row is for observed freq, 2nd row is for expected frequency
  int NumOfCol = this->data->vLocusInfo[iSnp].BothGenotypeCount.size();
  int totalGenotype = 0;
  double* contigency = new double[NumOfRow * NumOfCol];
  int idx = 0;
  for (genotype_iter = this->data->vLocusInfo[iSnp].BothGenotypeCount.begin();
       genotype_iter != this->data->vLocusInfo[iSnp].BothGenotypeCount.end();
       genotype_iter++) {
    double expectedFreq = GetExpectedGenotypeFreq(
        genotype_iter->first, this->data->vLocusInfo[iSnp].BothAlleleCount,
        AlleleType, this->vCoefficient, this->data->getSampleNum(),
        this->data->getNumOfChrSet());
    BOOST_ASSERT(idx < NumOfRow * NumOfCol);
    contigency[idx++] =
        genotype_iter
            ->second;  //(this->data->getCaseNum()*this->data->getNumOfChrSet());
    BOOST_ASSERT(idx < NumOfRow * NumOfCol);
    contigency[idx++] = expectedFreq;
    totalGenotype += genotype_iter->second;
  };

  // Fisher's exact test:
  BothChi = 0;
  for (int i = 0; i < NumOfCol * NumOfRow; i = i + 2) {
    if (contigency[i] < 1) contigency[i] *= totalGenotype;
    if (contigency[i + 1] < 1) contigency[i + 1] *= totalGenotype;
    if (contigency[i + 1] != 0)
    	BothChi += (contigency[i + 1] - contigency[i]) *
                 (contigency[i + 1] - contigency[i]) / contigency[i + 1];
  }
  boost::math::chi_squared dist(NumOfCol);
  BothPearsonP = boost::math::cdf(boost::math::complement(dist, BothChi));

  double expect = -1.0;
  double percnt = 100.0;
  double emin = 0;
  double pre = 0, prt = 0;
  int ws = 300000;
  try {
    fexact(&NumOfRow, &NumOfCol, contigency, &NumOfRow, &expect, &percnt, &emin,
           &prt, &pre, &ws);
    BothFisherP = pre;
  }
  catch (std::runtime_error&) {
    BothFisherP = -999;
  }
  // Pearson's ChiSquare test
  // PearsonChiSquareTest(contigency,NumOfRow,NumOfCol,CaseChi,CasePearsonP);
  delete[] contigency;
  contigency = 0;
}


std::string HWETest::reporthtml() {
  boost::shared_ptr<SHEsis::CreatHtmlTable> html(new SHEsis::CreatHtmlTable());
  html->createTable("HWETest");
  std::vector<std::string> data;
  data.push_back("SNP");
  if(this->data->vQuantitativeTrait.size()==0){
	  data.push_back("Chi<sup>2</sup> in case");
	  data.push_back("Pearson's p in case");
	  data.push_back("Fisher's p in case");
	  data.push_back("Chi<sup>2</sup> in ctrl");
	  data.push_back("Pearson's p in ctrl");
	  data.push_back("Fisher's p in ctrl");
  }
  data.push_back("Chi<sup>2</sup> in both");
  data.push_back("Pearson's p in both");
  data.push_back("Fisher's p in both");

  html->addHeadRow(data);
  for (int i = 0; i < this->vHWETestResult.size(); i++) {
    data.clear();
    data.push_back(this->data->vLocusName[i]);
    if(this->data->vQuantitativeTrait.size()==0){
		data.push_back(convert2string(this->vHWETestResult[i].CaseChiSquare));
		data.push_back(convert2string(this->vHWETestResult[i].CasePearsonP));
		data.push_back(convert2string(this->vHWETestResult[i].CaseFisherP));
		data.push_back(convert2string(this->vHWETestResult[i].ControlChiSquare));
		data.push_back(convert2string(this->vHWETestResult[i].ControlPearsonP));
		data.push_back(convert2string(this->vHWETestResult[i].ControlFisherP));
    };
    data.push_back(convert2string(this->vHWETestResult[i].BothChiSquare));
    data.push_back(convert2string(this->vHWETestResult[i].BothPearsonP));
    data.push_back(convert2string(this->vHWETestResult[i].BothFisherP));
    html->addDataRow(data);
  };
  return html->getTable();
}

void HWETest::printHWETestResults() {
  for (int i = 0; i < this->vHWETestResult.size(); i++) {
    std::cout << "\nLocus " << i
              << ":\nHWE test for case:\n(chi,pearsonp,fisherp)=("
              << this->vHWETestResult[i].CaseChiSquare << ","
              << this->vHWETestResult[i].CasePearsonP << ","
              << this->vHWETestResult[i].CaseFisherP << ")\n"
              << "HWE test for control:\n(chi,pearsonp,fisherp)=("
              << this->vHWETestResult[i].ControlChiSquare << ","
              << this->vHWETestResult[i].ControlPearsonP << ","
              << this->vHWETestResult[i].ControlFisherP << ")\n"
              << "HWE test for case and control:\n(chi,pearsonp,fisherp)=("
              << this->vHWETestResult[i].BothChiSquare << ","
              << this->vHWETestResult[i].BothPearsonP << ","
              << this->vHWETestResult[i].BothFisherP << ")\n";
  }
}

} /* namespace SHEsis */
