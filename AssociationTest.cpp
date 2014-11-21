/*
 * AssociationTest.cpp
 *
 *  Created on: Aug 8, 2014
 *      Author: ada
 */

#include "AssociationTest.h"
#include <boost/assert.hpp>
#include "utility.h"
#include "fisher.h"
#include <iostream>
#include <math.h>
#include <algorithm>
#include "CreatHtmlTable.h"
template std::string convert2string<double>(double v);
namespace SHEsis {

AssociationTest::AssociationTest(boost::shared_ptr<SHEsisData> mdata)
    : data(mdata),
      vAssocationTestResult(mdata->getSnpNum()),
      NumOfPermutation(-1),
      adjust(false) {};

AssociationTest::~AssociationTest() { vAssocationTestResult.clear(); }

std::string AssociationTest::reporttxt() {
  std::stringstream res;
  res.precision(3);
  res << "\n-------------------------------------------------------\n";
  res << "Single Locus Association Test (Binary phenotype)\n";
  res << "-------------------------------------------------------\n";
  for (int i = 0; i < this->vAssocationTestResult.size(); i++) {
    res << this->data->vLocusName[i] << "(Allele):\n";
    std::string detail = "\t";
    int count;
    boost::unordered_map<short, double>::iterator iter;
    count = 0;
    for (iter = this->data->vLocusInfo[i].BothAlleleCount.begin();
         iter != this->data->vLocusInfo[i].BothAlleleCount.end(); iter++) {
      detail += this->data->getallele(iter->first);
      if (count != this->data->vLocusInfo[i].BothAlleleCount.size() - 1)
        detail += "\t";
      count++;
    }
    detail += "\nCase\t";
    count = 0;
    for (iter = this->data->vLocusInfo[i].BothAlleleCount.begin();
         iter != this->data->vLocusInfo[i].BothAlleleCount.end(); iter++) {
      detail += convert2string(
                    this->data->vLocusInfo[i].CaseAlleleCount[iter->first]) +
                "(" +
                convert2string(
                    this->data->vLocusInfo[i].CaseAlleleCount[iter->first] /
                    (double)this->data->getCaseNum() /
                    this->data->getCaseCallrate(i) /
                    (double)this->data->getNumOfChrSet()) +
                ")";
      if (count != this->data->vLocusInfo[i].BothAlleleCount.size() - 1)
        detail += "\t";
      count++;
    }

    detail += "\nControl\t";
    count = 0;
    for (iter = this->data->vLocusInfo[i].BothAlleleCount.begin();
         iter != this->data->vLocusInfo[i].BothAlleleCount.end(); iter++) {
      detail += convert2string(
                    this->data->vLocusInfo[i].ControlAlleleCount[iter->first]) +
                "(" +
                convert2string(
                    this->data->vLocusInfo[i].ControlAlleleCount[iter->first] /
                    (double)this->data->getControlNum() /
                    this->data->getControlCallrate(i) /
                    (double)this->data->getNumOfChrSet()) +
                ")";
      if (count != this->data->vLocusInfo[i].BothAlleleCount.size() - 1)
        detail += "\t";
      count++;
    }
    res << detail << "\n";
    if (-999 != this->vAssocationTestResult[i].AlleleOddsRatio) {
      res << "Odds ratio=" << this->vAssocationTestResult[i].AlleleOddsRatio
          << ", 95% CI=["
          << this->vAssocationTestResult[i].AlleleOddsRatioLowLimit << "~"
          << this->vAssocationTestResult[i].AlleleOddsRatioUpLimit << "]\n";
    }
    res << "Call rate in cases is " << this->data->vLocusInfo[i].AlleleCallrate
        << "\n";
    res << "Chi2 is " << this->vAssocationTestResult[i].AlleleChiSquare << "\n";
    res << "Fisher's p value is "
        << convert2string(this->vAssocationTestResult[i].AlleleFisherP) << "\n";
    res << "Pearson's p value is "
        << convert2string(this->vAssocationTestResult[i].AllelePearsonP)
        << "\n";
    if (this->NumOfPermutation != -1) {
      res << "Permutation P value is "
          << convert2string(this->vAssocationTestResult[i].AllelePermutationP)
          << "\n";
    }
    if (this->adjust) {
      res << "Holm step-down adjusted p-value is "
          << convert2string(this->vAssocationTestResult[i].AlleleHolmP) << "\n";
      res << "Sidak single-step adjusted p-value is "
          << convert2string(this->vAssocationTestResult[i].AlleleSidakSSP)
          << "\n";
      res << "Sidak step-down adjusted p-value is "
          << convert2string(this->vAssocationTestResult[i].AlleleSidakSDP)
          << "\n";
      res << "Benjamini & Hochberg step-up FDR controled p value is "
          << convert2string(this->vAssocationTestResult[i].AlleleBHP) << "\n";
      res << "Benjamini & Yekutieli step-up FDR controled p value is "
          << convert2string(this->vAssocationTestResult[i].AlleleBYP) << "\n";
    }
    if (this->data->getNumOfChrSet() > 1) {
      res << "-------------------------------------------------------\n";
      res << this->data->vLocusName[i] << "(Genotype):\n";
      detail = "\t";
      boost::unordered_map<std::string, double>::iterator iter2;
      count = 0;
      for (iter2 = this->data->vLocusInfo[i].BothGenotypeCount.begin();
           iter2 != this->data->vLocusInfo[i].BothGenotypeCount.end();
           iter2++) {
        detail += this->data->getOriginGenotype(iter2->first);
        if (count != this->data->vLocusInfo[i].BothGenotypeCount.size() - 1)
          detail += "\t";
        count++;
      }
      detail += "\nCase\t";
      count = 0;
      for (iter2 = this->data->vLocusInfo[i].BothGenotypeCount.begin();
           iter2 != this->data->vLocusInfo[i].BothGenotypeCount.end();
           iter2++) {
        detail +=
            convert2string(
                this->data->vLocusInfo[i].CaseGenotypeCount[iter2->first]) +
            "(" +
            convert2string(
                this->data->vLocusInfo[i].CaseGenotypeCount[iter2->first] /
                (double)this->data->getCaseNum() /
                this->data->getCaseCallrate(i)) +
            ")";
        if (count != this->data->vLocusInfo[i].BothGenotypeCount.size() - 1)
          detail += "\t";
        count++;
      }

      detail += "\nControl\t";
      count = 0;
      for (iter2 = this->data->vLocusInfo[i].BothGenotypeCount.begin();
           iter2 != this->data->vLocusInfo[i].BothGenotypeCount.end();
           iter2++) {
        detail +=
            convert2string(
                this->data->vLocusInfo[i].ControlGenotypeCount[iter2->first]) +
            "(" +
            convert2string(
                this->data->vLocusInfo[i].ControlGenotypeCount[iter2->first] /
                (double)this->data->getControlNum() /
                this->data->getControlCallrate(i)) +
            ")";
        if (count != this->data->vLocusInfo[i].BothGenotypeCount.size() - 1)
          detail += "\t";
        count++;
      }
      res << detail << "\n";
      res << "Chi2 is " << this->vAssocationTestResult[i].GenotypeChiSquare
          << "\n";
      res << "Fisher's p is "
          << convert2string(this->vAssocationTestResult[i].GenotypeFisherP)
          << "\n";
      res << "Pearson's p is "
          << convert2string(this->vAssocationTestResult[i].GenoTypePearsonP)
          << "\n";
      if (this->NumOfPermutation != -1) {
        res << "Permutation P value is "
            << convert2string(
                   this->vAssocationTestResult[i].GenotypePermutationP) << "\n";
      }
      if (this->adjust) {
        res << "Holm step-down adjusted p-value is "
            << convert2string(this->vAssocationTestResult[i].GenotypeHolmP)
            << "\n";
        res << "Sidak single-step adjusted p-value is "
            << convert2string(this->vAssocationTestResult[i].GenotypeSidakSSP)
            << "\n";
        res << "Sidak step-down adjusted p-value is "
            << convert2string(this->vAssocationTestResult[i].GenotypeSidakSDP)
            << "\n";
        res << "Benjamini & Hochberg step-up FDR controled p value is "
            << convert2string(this->vAssocationTestResult[i].GenotypeBHP)
            << "\n";
        res << "Benjamini & Yekutieli step-up FDR controled p value is "
            << convert2string(this->vAssocationTestResult[i].GenotypeBYP)
            << "\n";
      }
    }
    res << "-------------------------------------------------------\n";
  }
  return res.str();
}

std::string AssociationTest::reporthtml() {
  std::stringstream res;
  res << "\n<h2> Single Locus Association Test (Binary phenotype): </h2>\n";
  res << "<h3>Alleles:</h3>\n";
  res << this->reporthtmlAllele();
  if (this->data->getNumOfChrSet() > 1) {
    res << "<h3>Genotypes:</h3>\n";
    res << this->reporthtmlGenotype();
  }
  return res.str();
}

std::string AssociationTest::reporthtmlAllele() {
  boost::shared_ptr<SHEsis::CreatHtmlTable> html(new SHEsis::CreatHtmlTable());
  html->createTable("Allele_Association");
  std::vector<std::string> data;
  data.push_back("SNP");
  data.push_back("Call rate");
  //  data.push_back("Allele Type");
  data.push_back("Chi<sup>2</sup>");
  data.push_back("Pearson's p");
  data.push_back("Fisher's p");
  if (this->NumOfPermutation != -1) data.push_back("Permutaion P");
  data.push_back("OR [95% CI]");
  if (this->adjust) {
    data.push_back("Holm");
    data.push_back("SidakSS");
    data.push_back("SidakSD");
    data.push_back("FDR_BH");
    data.push_back("FDR_BY");
  }
  data.push_back("Detail");
  html->addHeadRow(data);
  for (int i = 0; i < this->vAssocationTestResult.size(); i++) {
    data.clear();
    data.push_back(this->data->vLocusName[i]);
    data.push_back(convert2string(this->data->vLocusInfo[i].AlleleCallrate));
    //    std::string alleletype;
    std::string detail = "<pre>\n\t";
    int count;
    boost::unordered_map<short, double>::iterator iter;
    count = 0;
    for (iter = this->data->vLocusInfo[i].BothAlleleCount.begin();
         iter != this->data->vLocusInfo[i].BothAlleleCount.end(); iter++) {
      detail += this->data->getallele(iter->first);
      if (count != this->data->vLocusInfo[i].BothAlleleCount.size() - 1)
        detail += "\t";
      count++;
    }
    detail += "\nCase\t";
    count = 0;
    for (iter = this->data->vLocusInfo[i].BothAlleleCount.begin();
         iter != this->data->vLocusInfo[i].BothAlleleCount.end(); iter++) {
      detail += convert2string(
                    this->data->vLocusInfo[i].CaseAlleleCount[iter->first]) +
                "(" +
                convert2string(
                    this->data->vLocusInfo[i].CaseAlleleCount[iter->first] /
                    (double)this->data->getCaseNum() /
                    this->data->getCaseCallrate(i) /
                    (double)this->data->getNumOfChrSet()) +
                ")";
      if (count != this->data->vLocusInfo[i].BothAlleleCount.size() - 1)
        detail += "\t";
      count++;
    }

    detail += "\nControl\t";
    count = 0;
    for (iter = this->data->vLocusInfo[i].BothAlleleCount.begin();
         iter != this->data->vLocusInfo[i].BothAlleleCount.end(); iter++) {
      detail += convert2string(
                    this->data->vLocusInfo[i].ControlAlleleCount[iter->first]) +
                "(" +
                convert2string(
                    this->data->vLocusInfo[i].ControlAlleleCount[iter->first] /
                    (double)this->data->getControlNum() /
                    this->data->getControlCallrate(i) /
                    (double)this->data->getNumOfChrSet()) +
                ")";
      if (count != this->data->vLocusInfo[i].BothAlleleCount.size() - 1)
        detail += "\t";
      count++;
    }
    detail += "\n</pre>";

    //    int count = 0;
    //    for (iter = this->data->vLocusInfo[i].BothAlleleCount.begin();
    //         iter != this->data->vLocusInfo[i].BothAlleleCount.end(); iter++)
    // {
    ////      alleletype += this->data->getallele(iter->first);  //
    //
    //      detail += this->data->getallele(iter->first) + "(" +
    //                convert2string(
    //                    this->data->vLocusInfo[i].CaseAlleleCount[iter->first]
    // /
    //                    (double)this->data->getCaseNum() /
    //                    (double)this->data->getNumOfChrSet()) +
    //                "/" +
    //                convert2string(
    //                    this->data->vLocusInfo[i].ControlAlleleCount[iter->first]
    // /
    //                    (double)this->data->getControlNum() /
    //                    (double)this->data->getNumOfChrSet()) +
    //                ")";
    //      if (count != this->data->vLocusInfo[i].BothAlleleCount.size() - 1) {
    ////        alleletype += "/";
    //        detail += ", ";
    //      };
    //      count++;
    //    }
    //    data.push_back(alleletype);
    data.push_back(
        convert2string(this->vAssocationTestResult[i].AlleleChiSquare));
    data.push_back(
        convert2string(this->vAssocationTestResult[i].AllelePearsonP));
    data.push_back(
        convert2string(this->vAssocationTestResult[i].AlleleFisherP));
    if (this->NumOfPermutation != -1)
      data.push_back(
          convert2string(this->vAssocationTestResult[i].AllelePermutationP));

    std::string OR =
        convert2string(this->vAssocationTestResult[i].AlleleOddsRatio) + " [" +
        convert2string(this->vAssocationTestResult[i].AlleleOddsRatioLowLimit) +
        "~" +
        convert2string(this->vAssocationTestResult[i].AlleleOddsRatioUpLimit) +
        "]";
    data.push_back(OR);

    if (this->adjust) {
      data.push_back(
          convert2string(this->vAssocationTestResult[i].AlleleHolmP));
      data.push_back(
          convert2string(this->vAssocationTestResult[i].AlleleSidakSSP));
      data.push_back(
          convert2string(this->vAssocationTestResult[i].AlleleSidakSDP));
      data.push_back(convert2string(this->vAssocationTestResult[i].AlleleBHP));
      data.push_back(convert2string(this->vAssocationTestResult[i].AlleleBYP));
    }
    data.push_back(detail);
    html->addDataRow(data);
  };
  return html->getTable();
}

std::string AssociationTest::reporthtmlGenotype() {
  boost::shared_ptr<SHEsis::CreatHtmlTable> html(new SHEsis::CreatHtmlTable());
  html->createTable("Genotype_Association");
  std::vector<std::string> data;
  data.push_back("SNP");
  //  data.push_back("Call rate");
  data.push_back("Chi<sup>2</sup>");
  data.push_back("Pearson's p");
  data.push_back("Fisher's p");
  if (this->NumOfPermutation != -1) data.push_back("Permutaion P");
  if (this->adjust) {
    data.push_back("Holm");
    data.push_back("SidakSS");
    data.push_back("SidakSD");
    data.push_back("FDR_BH");
    data.push_back("FDR_BY");
  }
  data.push_back("Detail");
  html->addHeadRow(data);
  for (int i = 0; i < this->vAssocationTestResult.size(); i++) {
    data.clear();
    data.push_back(this->data->vLocusName[i]);
    data.push_back(
        convert2string(this->vAssocationTestResult[i].GenotypeChiSquare));
    data.push_back(
        convert2string(this->vAssocationTestResult[i].GenoTypePearsonP));
    data.push_back(
        convert2string(this->vAssocationTestResult[i].GenotypeFisherP));
    if (this->NumOfPermutation != -1)
      data.push_back(
          convert2string(this->vAssocationTestResult[i].GenotypePermutationP));
    if (this->adjust) {
      data.push_back(
          convert2string(this->vAssocationTestResult[i].GenotypeHolmP));
      data.push_back(
          convert2string(this->vAssocationTestResult[i].GenotypeSidakSSP));
      data.push_back(
          convert2string(this->vAssocationTestResult[i].GenotypeSidakSDP));
      data.push_back(
          convert2string(this->vAssocationTestResult[i].GenotypeBHP));
      data.push_back(
          convert2string(this->vAssocationTestResult[i].GenotypeBYP));
    }

    std::string detail = "<pre>\n\t";
    boost::unordered_map<std::string, double>::iterator iter;
    int count = 0;
    for (iter = this->data->vLocusInfo[i].BothGenotypeCount.begin();
         iter != this->data->vLocusInfo[i].BothGenotypeCount.end(); iter++) {
      detail += this->data->getOriginGenotype(iter->first);
      if (count != this->data->vLocusInfo[i].BothGenotypeCount.size() - 1)
        detail += "\t";
      count++;
    }
    detail += "\nCase\t";
    count = 0;
    for (iter = this->data->vLocusInfo[i].BothGenotypeCount.begin();
         iter != this->data->vLocusInfo[i].BothGenotypeCount.end(); iter++) {
      detail += convert2string(
                    this->data->vLocusInfo[i].CaseGenotypeCount[iter->first]) +
                "(" +
                convert2string(
                    this->data->vLocusInfo[i].CaseGenotypeCount[iter->first] /
                    (double)this->data->getCaseNum() /
                    this->data->getCaseCallrate(i)) +
                ")";
      if (count != this->data->vLocusInfo[i].BothGenotypeCount.size() - 1)
        detail += "\t";
      count++;
    }

    detail += "\nControl\t";
    count = 0;
    for (iter = this->data->vLocusInfo[i].BothGenotypeCount.begin();
         iter != this->data->vLocusInfo[i].BothGenotypeCount.end(); iter++) {
      detail +=
          convert2string(
              this->data->vLocusInfo[i].ControlGenotypeCount[iter->first]) +
          "(" +
          convert2string(
              this->data->vLocusInfo[i].ControlGenotypeCount[iter->first] /
              (double)this->data->getControlNum() /
              this->data->getControlCallrate(i)) +
          ")";
      if (count != this->data->vLocusInfo[i].BothGenotypeCount.size() - 1)
        detail += "\t";
      count++;
    }
    detail += "\n</pre>";
    //    for (iter = this->data->vLocusInfo[i].BothGenotypeCount.begin();
    //         iter != this->data->vLocusInfo[i].BothGenotypeCount.end();
    // iter++) {
    //      detail +=
    //          this->data->getOriginGenotype(iter->first) + "(" +
    //          convert2string(
    //              this->data->vLocusInfo[i].CaseGenotypeCount[iter->first] /
    //              (double)this->data->getCaseNum()) +
    //          "/" +
    //          convert2string(
    //              this->data->vLocusInfo[i].ControlGenotypeCount[iter->first]
    // /
    //              (double)this->data->getControlNum()) +
    //          ")";
    //      if (count != this->data->vLocusInfo[i].BothGenotypeCount.size() - 1)
    // {
    //        detail += ", ";
    //      };
    //      count++;
    //    }

    data.push_back(detail);
    html->addDataRow(data);
  };
  return html->getTable();
}

void AssociationTest::association() {
  if (this->data->vLocusInfo[0].BothAlleleCount.size() == 0 &&
      this->data->vQuantitativeTrait.size() > 0)
    this->data->statCount();
  else if (this->data->vLocusInfo[0].BothAlleleCount.size() == 0 &&
           this->data->vQuantitativeTrait.size() == 0)
    this->data->statCount(this->data->vLabel);
  this->AssociationTestForAllSnpsAllele(false);
  if (this->data->getNumOfChrSet() > 1)
    this->AssociationTestForAllSnpsGenotype(false);
}

void AssociationTest::AssociationTestForAllSnpsAllele(
    bool bPermutating = false) {
  for (int i = 0; i < this->vAssocationTestResult.size(); i++) {
    this->SingleSnpTestAllele(
        i, this->vAssocationTestResult[i].AlleleFisherP,
        this->vAssocationTestResult[i].AllelePearsonP,
        this->vAssocationTestResult[i].AlleleChiSquare,
        this->vAssocationTestResult[i].AlleleOddsRatio,
        this->vAssocationTestResult[i].AlleleOddsRatioLowLimit,
        this->vAssocationTestResult[i].AlleleOddsRatioUpLimit, bPermutating);
  }
  // adjust p
  if (bPermutating) return;
  if (this->adjust) {
    std::vector<MultiComp> originp;
    std::vector<double> adjusted;
    for (int i = 0; i < this->vAssocationTestResult.size(); i++) {
      MultiComp val;
      val.p = this->vAssocationTestResult[i].AllelePearsonP;
      val.idx = i;
      originp.push_back(val);
    }

    std::sort(originp.begin(), originp.end());
    HolmCorrection(originp, adjusted);
    for (int i = 0; i < this->vAssocationTestResult.size(); i++)
      this->vAssocationTestResult[i].AlleleHolmP = adjusted[i];

    SidakSDCorrection(originp, adjusted);
    for (int i = 0; i < this->vAssocationTestResult.size(); i++)
      this->vAssocationTestResult[i].AlleleSidakSDP = adjusted[i];

    SidakSSCorrection(originp, adjusted);
    for (int i = 0; i < this->vAssocationTestResult.size(); i++)
      this->vAssocationTestResult[i].AlleleSidakSSP = adjusted[i];

    BHCorrection(originp, adjusted);
    for (int i = 0; i < this->vAssocationTestResult.size(); i++)
      this->vAssocationTestResult[i].AlleleBHP = adjusted[i];

    BYCorrection(originp, adjusted);
    for (int i = 0; i < this->vAssocationTestResult.size(); i++)
      this->vAssocationTestResult[i].AlleleBYP = adjusted[i];
  }
}

void AssociationTest::AssociationTestForAllSnpsGenotype(
    bool bPermutating = false) {
  for (int i = 0; i < this->vAssocationTestResult.size(); i++) {
    this->SingleSnpTestGenotype(
        i, this->vAssocationTestResult[i].GenotypeFisherP,
        this->vAssocationTestResult[i].GenoTypePearsonP,
        this->vAssocationTestResult[i].GenotypeChiSquare, bPermutating);
  }
  // adjust p
  if (bPermutating) return;
  if (this->adjust) {
    // allele
    std::vector<MultiComp> originp;
    std::vector<double> adjusted;
    for (int i = 0; i < this->vAssocationTestResult.size(); i++) {
      MultiComp val;
      val.p = this->vAssocationTestResult[i].GenoTypePearsonP;
      val.idx = i;
      originp.push_back(val);
    };
    std::sort(originp.begin(), originp.end());
    HolmCorrection(originp, adjusted);
    for (int i = 0; i < this->vAssocationTestResult.size(); i++)
      this->vAssocationTestResult[i].GenotypeHolmP = adjusted[i];

    SidakSDCorrection(originp, adjusted);
    for (int i = 0; i < this->vAssocationTestResult.size(); i++)
      this->vAssocationTestResult[i].GenotypeSidakSDP = adjusted[i];

    SidakSSCorrection(originp, adjusted);
    for (int i = 0; i < this->vAssocationTestResult.size(); i++)
      this->vAssocationTestResult[i].GenotypeSidakSSP = adjusted[i];

    BHCorrection(originp, adjusted);
    for (int i = 0; i < this->vAssocationTestResult.size(); i++)
      this->vAssocationTestResult[i].GenotypeBHP = adjusted[i];

    BYCorrection(originp, adjusted);
    for (int i = 0; i < this->vAssocationTestResult.size(); i++)
      this->vAssocationTestResult[i].GenotypeBYP = adjusted[i];
  }
}

void getTheSmallestP(std::vector<LocusAssiciationTestResult> res,
                     double& allelep, double& genop) {
  allelep = 1;
  genop = 1;
  for (int i = 0; i < res.size(); i++) {
    if (res[i].AllelePearsonP < allelep) {
      allelep = res[i].AllelePearsonP;
    };
    if (res[i].GenoTypePearsonP < genop) {
      genop = res[i].GenoTypePearsonP;
    }
  }
}

// int getRank(double p, std::vector<double> v) {
//  for (int i = 0; i < v.size() - 1; i++) {
//    if (p >= v[i] && p <= v[i + 1]) return i;
//  }
//  return v.size();
//}

void AssociationTest::permutation() {
  this->vPermutateLabel = this->data->vLabel;
  for (int i = 0; i < this->NumOfPermutation; i++) {
    if (i % 200 == 0) {
      printf("\rPermutating...%d%%",
             (int)(100 * (double)i / (double)this->NumOfPermutation));
      fflush(stdout);
    }
    std::random_shuffle(this->vPermutateLabel.begin(),
                        this->vPermutateLabel.end());
    this->data->statCount(this->vPermutateLabel);
    this->AssociationTestForAllSnpsAllele(true);
    if (this->data->getNumOfChrSet() > 1)
      this->AssociationTestForAllSnpsGenotype(true);
    double ap, gp;
    getTheSmallestP(this->vAssocationTestResult, ap, gp);
    this->PermutationPAllele.push_back(ap);
    this->PermutationPGenotype.push_back(gp);
  };
  printf("\rPermutating...%d%%\n", 100);
  fflush(stdout);
  std::sort(this->PermutationPAllele.begin(), this->PermutationPAllele.end());
  std::sort(this->PermutationPGenotype.begin(),
            this->PermutationPGenotype.end());
  //  for(int i=0;i<this->PermutationPAllele.size();i++){
  //	  std::cout<<this->PermutationPAllele[i]<<",";
  //  }
  this->data->statCount(data->vLabel);
  std::cout << "\n";
  this->AssociationTestForAllSnpsAllele(false);
  if (this->data->getNumOfChrSet() > 1)
    this->AssociationTestForAllSnpsGenotype(false);
  for (int i = 0; i < this->vAssocationTestResult.size(); i++) {
    this->vAssocationTestResult[i].AllelePermutationP =
        this->vAssocationTestResult[i].AllelePearsonP > 0
            ? (double)getRank(this->vAssocationTestResult[i].AllelePearsonP,
                              this->PermutationPAllele) /
                  (double)this->NumOfPermutation
            : -999;
    if (this->data->getNumOfChrSet() > 1)
      this->vAssocationTestResult[i].GenotypePermutationP =
          this->vAssocationTestResult[i].GenoTypePearsonP > 0
              ? (double)getRank(this->vAssocationTestResult[i].GenoTypePearsonP,
                                this->PermutationPGenotype) /
                    (double)this->NumOfPermutation
              : -999;
    //    std::cout<<"allep,genop="<<this->vAssocationTestResult[i].AllelePearsonP<<","<<this->vAssocationTestResult[i].GenoTypePearsonP<<"\n";
  };
}

void AssociationTest::printAssociationTestResults() {
  for (int i = 0; i < this->vAssocationTestResult.size(); i++) {
    std::cout << "\nLocus " << i << "\nAllele association test:\n"
              << "(fisher's p, pearson's p, chi suqare, permutation p,or, "
                 "orLow, orUp)="
              << "(" << vAssocationTestResult[i].AlleleFisherP << ","
              << vAssocationTestResult[i].AllelePearsonP << ","
              << vAssocationTestResult[i].AlleleChiSquare << ","
              << vAssocationTestResult[i].AllelePermutationP << ","
              << vAssocationTestResult[i].AlleleOddsRatio << ","
              << vAssocationTestResult[i].AlleleOddsRatioLowLimit << ","
              << vAssocationTestResult[i].AlleleOddsRatioUpLimit << ")\n";
    std::cout << "Genotype association test:\n"
              << "(fisher's p, pearson's p, chi square,permutation p)=("
              << vAssocationTestResult[i].GenotypeFisherP << ","
              << vAssocationTestResult[i].GenoTypePearsonP << ","
              << vAssocationTestResult[i].GenotypeChiSquare << ","
              << vAssocationTestResult[i].GenotypePermutationP << ")\n";
  }
}

void AssociationTest::SingleSnpTestAllele(int iSnp, double& FisherP,
                                          double& PearsonP, double& ChiSquare,
                                          double& oddsRatio, double& ORLowLimit,
                                          double& ORUpLimit,
                                          bool bPermutating = false) {

  BOOST_ASSERT(this->data->vLocusInfo[iSnp].CaseAlleleCount.size() ==
               this->data->vLocusInfo[iSnp].ControlAlleleCount.size());

  int NumOfAlleleType = this->data->vLocusInfo[iSnp].ControlAlleleCount.size();
  int NumOfPhenotype = 2;
  double* contigency = new double[NumOfPhenotype * NumOfAlleleType];

  boost::unordered_map<short, double>::iterator map_it;
  int idx = 0;
  for (map_it = this->data->vLocusInfo[iSnp].CaseAlleleCount.begin();
       map_it != this->data->vLocusInfo[iSnp].CaseAlleleCount.end(); map_it++) {
    BOOST_ASSERT(
        this->data->vLocusInfo[iSnp].ControlAlleleCount.end() !=
        this->data->vLocusInfo[iSnp].ControlAlleleCount.find(map_it->first));
    contigency[idx++] =
        this->data->vLocusInfo[iSnp].ControlAlleleCount[map_it->first];
    contigency[idx++] = map_it->second;
  };

  // Pearson's ChiSquare test
  PearsonChiSquareTest(contigency, NumOfPhenotype, NumOfAlleleType, ChiSquare,
                       PearsonP);

  if (!bPermutating) {
    // Fisher's exact test:
    double expect = -1.0;
    double percnt = 100.0;
    double emin = 0;
    double pre = 0, prt = 0;
    int ws = 300000;
    try {
      fexact(&NumOfPhenotype, &NumOfAlleleType, contigency, &NumOfPhenotype,
             &expect, &percnt, &emin, &prt, &pre, &ws);
      FisherP = pre;
    }
    catch (std::runtime_error&) {
      FisherP = -999;
    }

    // get odds ratio
    if ((2 == NumOfAlleleType) && (0 != contigency[3] && 0 != contigency[0] &&
                                   0 != contigency[1] && 0 != contigency[2])) {
      oddsRatio =
          (contigency[1] * contigency[2] / (contigency[0] * contigency[3]));
      // output by maf
      map_it = this->data->vLocusInfo[iSnp].ControlAlleleCount.begin();
      double ctrlfreq = map_it->second / (double)this->data->getControlNum() /
                        (double)this->data->getNumOfChrSet();
      oddsRatio = ctrlfreq > 0.5 ? 1 / oddsRatio : oddsRatio;
      double v = 1 / contigency[1] + 1 / contigency[2] + 1 / contigency[3] +
                 1 / contigency[0];

      ORLowLimit = oddsRatio * exp(-1.96 * sqrt(v));
      ORUpLimit = oddsRatio * exp(1.96 * sqrt(v));
    } else {
      oddsRatio = -999;
      ORLowLimit = -999;
      ORUpLimit = -999;
    }
  }
  delete[] contigency;
  contigency = 0;
}

void AssociationTest::SingleSnpTestGenotype(int iSnp, double& FisherP,
                                            double& PearsonP, double& ChiSquare,
                                            bool bPermutating = false) {

  BOOST_ASSERT(this->data->vLocusInfo[iSnp].CaseGenotypeCount.size() ==
               this->data->vLocusInfo[iSnp].ControlGenotypeCount.size());

  int NumOfAlleleType =
      this->data->vLocusInfo[iSnp].ControlGenotypeCount.size();
  int NumOfPhenotype = 2;
  double* contigency = new double[NumOfPhenotype * NumOfAlleleType];

  boost::unordered_map<std::string, double>::iterator map_it;
  int idx = 0;
  for (map_it = this->data->vLocusInfo[iSnp].CaseGenotypeCount.begin();
       map_it != this->data->vLocusInfo[iSnp].CaseGenotypeCount.end();
       map_it++) {
    BOOST_ASSERT(
        this->data->vLocusInfo[iSnp].ControlGenotypeCount.end() !=
        this->data->vLocusInfo[iSnp].ControlGenotypeCount.find(map_it->first));
    contigency[idx++] =
        this->data->vLocusInfo[iSnp].ControlGenotypeCount[map_it->first];
    contigency[idx++] = map_it->second;
  };
  if (!bPermutating) {
    // Fisher's exact test:
    double expect = -1.0;
    double percnt = 100.0;
    double emin = 0;
    double pre = 0, prt = 0;
    int ws = 300000;
    try {
      fexact(&NumOfPhenotype, &NumOfAlleleType, contigency, &NumOfPhenotype,
             &expect, &percnt, &emin, &prt, &pre, &ws);
      FisherP = pre;
    }
    catch (std::runtime_error&) {
      FisherP = -999;
    }
  }
  // Pearson's ChiSquare test
  PearsonChiSquareTest(contigency, NumOfPhenotype, NumOfAlleleType, ChiSquare,
                       PearsonP);
  delete[] contigency;
  contigency = 0;
}

} /* namespace SHEsis */
