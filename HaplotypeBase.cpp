/*
 * HaplotypeBase.cpp
 *
 *  Created on: Sep 20, 2014
 *      Author: ionadmin
 */
#include "HaplotypeBase.h"
#include "CreatHtmlTable.h"
#include <math.h>
#include <boost/math/distributions/students_t.hpp>
template std::string convert2string<double>(double v);
namespace SHEsis {
int HaplotypeBase::getHaploCount(int sample, int hapIdx) {
  int count = 0;
  for (int p = 0; p < this->data->getNumOfChrSet(); p++) {
    if (this->Results.genotypes[sample][p] == hapIdx) count++;
  }
  return count;
}
bool HaplotypeBase::IsHaploMissing(int sample) {
  for (int p = 0; p < this->data->getNumOfChrSet(); p++) {
    if (this->Results.genotypes[sample][p] == -1) return true;
  }
  return false;
}
singHapQTLRes HaplotypeBase::SingleHaploAssociationTestQTL(int hapIdx) {
  singHapQTLRes res;
  //    double freq = ((double)this->Results.BothCount[hapIdx])/
  //                  ((double)this->data->getSampleNum() *
  //                   (double)this->data->getNumOfChrSet());
  double freq = ((double)this->Results.BothCount[hapIdx]) /
                ((double)(this->NonmissingCase + this->NonmissingCtrl) *
                 (double)this->data->getNumOfChrSet());
  //    std::cout<<"freq="<<freq<<", nonmissing
  // case="<<this->NonmissingCase<<",nonmissing
  // ctrl="<<this->NonmissingCtrl<<"\n";
  if (freq < this->freqthreshold) return res;

  double g_mean = 0, g_var = 0;
  double qt_mean = 0, qt_var = 0;
  double qt_g_covar = 0;
  int ValidSampleNum = 0;
  for (int iSample = 0; iSample < this->data->getSampleNum(); iSample++) {
    qt_mean += this->data->vQuantitativeTrait[iSample];
    g_mean += this->getHaploCount(iSample, hapIdx);
    ValidSampleNum++;
  }
  qt_mean /= (double)ValidSampleNum;
  g_mean /= (double)ValidSampleNum;

  for (int iSample = 0; iSample < this->data->getSampleNum(); iSample++) {
    if (0 == this->data->vQuantitativeTrait[iSample]) continue;
    if (this->IsHaploMissing(iSample)) continue;
    qt_var += (this->data->vQuantitativeTrait[iSample] - qt_mean) *
              (this->data->vQuantitativeTrait[iSample] - qt_mean);
    double g = this->getHaploCount(iSample, hapIdx);
    g_var += (g - g_mean) * (g - g_mean);
    qt_g_covar +=
        (this->data->vQuantitativeTrait[iSample] - qt_mean) * (g - g_mean);
  }

  qt_var /= (double)ValidSampleNum - 1;
  g_var /= (double)ValidSampleNum - 1;
  qt_g_covar /= (double)ValidSampleNum - 1;
  res.beta = qt_g_covar / g_var;
  res.se = sqrt((qt_var / g_var - (qt_g_covar * qt_g_covar) / (g_var * g_var)) /
                ((double)ValidSampleNum - 2));
  res.T = res.beta / res.se;
  if (res.T < 0.0000001) res.T = 0;
  if (res.T > 99999) res.T = 99999;
  if (ValidSampleNum - 2 < 1)
    res.p = -999;
  else {
    boost::math::students_t dist(ValidSampleNum - 2);
    res.p = 2 * boost::math::cdf(boost::math::complement(dist, fabs(res.T)));
  }
  res.R2 = (qt_g_covar * qt_g_covar) / (qt_var * g_var);
  res.ValidSampleNum = ValidSampleNum;
  return res;
}

void HaplotypeBase::AssociationTestQTL() {
  for (int hapIdx = 0; hapIdx < this->Results.haplotypes.size(); hapIdx++) {
    singHapQTLRes res = SingleHaploAssociationTestQTL(hapIdx);
    this->Results.singleHapQTL.push_back(res);
  }
  if (this->adjust) {
    std::vector<MultiComp> originp;
    std::vector<double> adjusted;
    for (int i = 0; i < this->Results.singleHapQTL.size(); i++) {
      MultiComp val;
      val.p = this->Results.singleHapQTL[i].p;
      val.idx = i;
      originp.push_back(val);
    }
    std::sort(originp.begin(), originp.end());
    HolmCorrection(originp, adjusted);
    for (int i = 0; i < this->Results.singleHapQTL.size(); i++)
      this->Results.singleHapQTL[i].HolmP = adjusted[i];

    SidakSDCorrection(originp, adjusted);
    for (int i = 0; i < this->Results.singleHapQTL.size(); i++)
      this->Results.singleHapQTL[i].SidakSDP = adjusted[i];

    SidakSSCorrection(originp, adjusted);
    for (int i = 0; i < this->Results.singleHapQTL.size(); i++)
      this->Results.singleHapQTL[i].SidakSSP = adjusted[i];

    BHCorrection(originp, adjusted);
    for (int i = 0; i < this->Results.singleHapQTL.size(); i++)
      this->Results.singleHapQTL[i].BHP = adjusted[i];

    BYCorrection(originp, adjusted);
    for (int i = 0; i < this->Results.singleHapQTL.size(); i++)
      this->Results.singleHapQTL[i].BYP = adjusted[i];
  };
}

void HaplotypeBase::AssociationTest() {

  if (this->data->vQuantitativeTrait.size() == 0) {
    this->AssociationTestBinary();
  } else {
    this->AssociationTestQTL();
  }
}

void HaplotypeBase::AssociationTestBinary() {

  int haploNum = this->Results.haplotypes.size();
  this->Results.singleHap.resize(haploNum);
  int nrow = 2;
  double expect = -1.0;
  double percnt = 100.0;
  double emin = 0;
  double pre = 0, prt = 0;
  int ws = 300000;

  int validHap = 0;
  double* contigency = new double[4];
  // single hapotype statistics
  for (int i = 0; i < haploNum; i++) {
    double freq = ((double)this->Results.ControlCount[i] +
                   (double)this->Results.CaseCount[i]) /
                  ((double)(this->NonmissingCase + this->NonmissingCtrl) *
                   (double)this->data->getNumOfChrSet());
    if (freq < this->freqthreshold) continue;
    validHap++;
    contigency[0] = this->Results.ControlCount[i];
    contigency[1] = this->Results.CaseCount[i];
    contigency[2] = this->data->getControlNum() * this->data->getNumOfChrSet() -
                    contigency[0];
    contigency[3] =
        this->data->getCaseNum() * this->data->getNumOfChrSet() - contigency[1];
    try {
      fexact(&nrow, &nrow, contigency, &nrow, &expect, &percnt, &emin, &prt,
             &pre, &ws);
      this->Results.singleHap[i].fisherp = pre;
    }
    catch (std::runtime_error&) {
      this->Results.singleHap[i].fisherp = -999;
    }
    PearsonChiSquareTest(contigency, nrow, nrow,
                         this->Results.singleHap[i].chisquare,
                         this->Results.singleHap[i].pearsonp);
    if ((0 != contigency[3] && 0 != contigency[0] && 0 != contigency[1] &&
         0 != contigency[2])) {
      this->Results.singleHap[i].OR =
          (contigency[1] * contigency[2] / (contigency[0] * contigency[3]));
      double v = 1 / contigency[1] + 1 / contigency[2] + 1 / contigency[3] +
                 1 / contigency[0];
      this->Results.singleHap[i].orlow =
          this->Results.singleHap[i].OR * exp(-1.96 * sqrt(v));
      this->Results.singleHap[i].orUp =
          this->Results.singleHap[i].OR * exp(1.96 * sqrt(v));
    } else {
      this->Results.singleHap[i].OR = -999;
      this->Results.singleHap[i].orlow = -999;
      this->Results.singleHap[i].orUp = -999;
    }
  }

  delete[] contigency;
  if (validHap == 0) return;
  contigency = new double[2 * validHap];
  int idx = 0;
  int HapCount = 0;
  for (int i = 0; i < haploNum; i++) {
    double freq = ((double)this->Results.ControlCount[i] +
                   (double)this->Results.CaseCount[i]) /
                  ((double)(this->NonmissingCase + this->NonmissingCtrl) *
                   (double)this->data->getNumOfChrSet());
    if (freq < this->freqthreshold) continue;
    HapCount++;
    contigency[idx++] = this->Results.ControlCount[i];
    contigency[idx++] = this->Results.CaseCount[i];
  };

  try {
    fexact(&nrow, &validHap, contigency, &nrow, &expect, &percnt, &emin, &prt,
           &pre, &ws);
    this->Results.FisherP = pre;
  }
  catch (std::runtime_error&) {
    this->Results.FisherP = -999;
  }
  // Pearson's ChiSquare test
  PearsonChiSquareTest(contigency, nrow, validHap, this->Results.ChiSquare,
                       this->Results.PearsonP);
  delete[] contigency;

  if (this->adjust) {
    std::vector<MultiComp> originp;
    std::vector<double> adjusted;
    for (int i = 0; i < this->Results.singleHap.size(); i++) {
      MultiComp val;
      val.p = this->Results.singleHap[i].pearsonp;
      val.idx = i;
      originp.push_back(val);
    };

    std::sort(originp.begin(), originp.end());
    HolmCorrection(originp, adjusted);
    for (int i = 0; i < this->Results.singleHap.size(); i++)
      this->Results.singleHap[i].HolmP = adjusted[i];

    SidakSDCorrection(originp, adjusted);
    for (int i = 0; i < this->Results.singleHap.size(); i++)
      this->Results.singleHap[i].SidakSDP = adjusted[i];

    SidakSSCorrection(originp, adjusted);
    for (int i = 0; i < this->Results.singleHap.size(); i++)
      this->Results.singleHap[i].SidakSSP = adjusted[i];

    BHCorrection(originp, adjusted);
    for (int i = 0; i < this->Results.singleHap.size(); i++)
      this->Results.singleHap[i].BHP = adjusted[i];

    BYCorrection(originp, adjusted);
    for (int i = 0; i < this->Results.singleHap.size(); i++)
      this->Results.singleHap[i].BYP = adjusted[i];
  }
}

std::string HaplotypeBase::reporthtml() {
  std::stringstream ss;
  ss << "\n<h2> Haplotype Analysis: </h2>\n";
  ss << "<p>Haplotypes with frequency <" << this->freqthreshold
     << " are ignored.<br>";
  ss << "Loci chosen for haplotype analysis: ";
  for (int i = 0; i < this->SnpIdx.size() - 1; i++) {
    ss << this->data->vLocusName[SnpIdx[i]] << ", ";
  }
  ss << this->data->vLocusName[SnpIdx[this->SnpIdx.size() - 1]] << "</p>\n";
  if (this->data->vQuantitativeTrait.size() == 0) {
    ss << this->reporthtmltableBinary();
    ss << "<p><b>Global result:</b><br>Total control="
       << this->data->getControlNum()
       << ", total case=" << this->data->getCaseNum() << ".<br>";
    ss << "Global Chi<sup>2</sup> is "
       << convert2string(this->Results.ChiSquare) << ", ";
    ss << "Fisher's p is " << convert2string(this->Results.FisherP) << ", ";
    ss << "Pearson's p is " << convert2string(this->Results.PearsonP)
       << ".</p>\n";
  } else {
    ss << this->reporthtmltableQTL();
  }
  return ss.str();
}

std::string HaplotypeBase::reporttxt() {
  std::stringstream res;
  res.precision(3);
  res << "\n-------------------------------------------------------\n";
  res << "Haplotype Analysis\n";
  res << "-------------------------------------------------------\n";
  res << "Haplotypes with frequency <" << this->freqthreshold
      << " are ignored.\n";
  res << "Loci chosen for haplotype analysis: ";
  for (int i = 0; i < this->SnpIdx.size() - 1; i++) {
    res << this->data->vLocusName[SnpIdx[i]] << ", ";
  };
  res << this->data->vLocusName[SnpIdx[this->SnpIdx.size() - 1]] << "\n";
  res << "-------------\n";
  if (this->data->vQuantitativeTrait.size() == 0) {
    res << this->reporttxttableBinary();
    res << "-------------\n";
    res << "Global result:\nTotal control=" << this->data->getControlNum()
        << ", total case=" << this->data->getCaseNum() << ".\n";
    res << "Global Chi2 is " << convert2string(this->Results.ChiSquare) << ", ";
    res << "Fisher's p is " << convert2string(this->Results.FisherP) << ", ";
    res << "Pearson's p is " << convert2string(this->Results.PearsonP) << ".\n";
  } else {
    res << this->reporttxttableQTL();
  }
  res << "-------------------------------------------------------\n";
  return res.str();
}

std::string HaplotypeBase::reporttxttableBinary() {
  std::stringstream res;
  res.precision(3);
  res << "Haplotype\tCase(freq)\tControl(freq)\tChi2\t\tFisher's p\tPearson's "
         "p\tOR [95% CI]";
  if (this->adjust) {
    res << "\tHolm\tSidakSS\tSidakSD\tFDR_BH\tFDR_BY";
  }
  res << "\n";
  for (int i = 0; i < this->Results.singleHap.size(); i++) {
    if (-999 == this->Results.singleHap[i].pearsonp) continue;
    for (int j = 0; j < this->SnpIdx.size(); j++) {
      res << this->data->getallele(this->Results.haplotypes[i][j]);
    }
    res << "\t\t";
    double freq =
        (double)this->Results.CaseCount[i] /
        //                ((double)this->data->getCaseNum() *
        ((double)this->NonmissingCase * (double)this->data->getNumOfChrSet());
    res << this->Results.CaseCount[i] << "(" << freq << ")\t\t";
    freq =
        (double)this->Results.ControlCount[i] /
        //((double)this->data->getControlNum() *
        ((double)this->NonmissingCtrl * (double)this->data->getNumOfChrSet());
    res << this->Results.ControlCount[i] << "(" << freq << ")\t\t";
    res << this->Results.singleHap[i].chisquare << "\t\t";
    res << this->Results.singleHap[i].pearsonp << "\t\t";
    res << convert2string(this->Results.singleHap[i].fisherp) << "\t\t";
    res << convert2string(this->Results.singleHap[i].OR) << "["
        << convert2string(this->Results.singleHap[i].orlow) << "~"
        << convert2string(this->Results.singleHap[i].orUp) << "]";
    if (this->adjust) {
      res << "\t" << this->Results.singleHap[i].HolmP << "\t"
          << this->Results.singleHap[i].SidakSSP << "\t"
          << this->Results.singleHap[i].SidakSDP << "\t"
          << this->Results.singleHap[i].BHP << "\t"
          << this->Results.singleHap[i].BYP;
    }
    res << "\n";
  };
  return res.str();
}

std::string HaplotypeBase::reporthtmltableBinary() {
  boost::shared_ptr<SHEsis::CreatHtmlTable> html(new SHEsis::CreatHtmlTable());
  html->createTable("Haplotype_Analysis");
  std::vector<std::string> data;
  data.push_back("Haplotype");
  data.push_back("Case(freq)");
  data.push_back("Control(freq)");
  data.push_back("Chi<sup>2</sup>");
  data.push_back("Fisher's p");
  data.push_back("Pearson's p");
  data.push_back("OR [95% CI]");
  if (this->adjust) {
    data.push_back("Holm");
    data.push_back("SidakSS");
    data.push_back("SidakSD");
    data.push_back("FDR_BH");
    data.push_back("FDR_BY");
  }
  html->addHeadRow(data);
  for (int i = 0; i < this->Results.singleHap.size(); i++) {
    if (-999 == this->Results.singleHap[i].pearsonp) continue;
    data.clear();
    std::stringstream hap;
    for (int j = 0; j < this->SnpIdx.size(); j++) {
      hap << this->data->getallele(this->Results.haplotypes[i][j]);
    }
    data.push_back(hap.str());
    std::string casefreq, controlfreq;
    casefreq +=
        convert2string((double)this->Results.CaseCount[i]) + "(" +
        convert2string(
            (double)this->Results.CaseCount[i] /
            //                               ((double)this->data->getCaseNum() *
            ((double)this->NonmissingCase *
             (double)this->data->getNumOfChrSet())) +
        ")";
    data.push_back(casefreq);
    controlfreq +=
        convert2string((double)this->Results.ControlCount[i]) + "(" +
        convert2string(
            (double)this->Results.ControlCount[i] /
            //                                  ((double)this->data->getControlNum()
            // *
            ((double)this->NonmissingCtrl *
             (double)this->data->getNumOfChrSet())) +
        ")";
    data.push_back(controlfreq);
    data.push_back(convert2string(this->Results.singleHap[i].chisquare));
    data.push_back(convert2string(this->Results.singleHap[i].fisherp));
    data.push_back(convert2string(this->Results.singleHap[i].pearsonp));

    std::string OR;
    OR += convert2string(this->Results.singleHap[i].OR) + " [" +
          convert2string(this->Results.singleHap[i].orlow) + "~" +
          convert2string(this->Results.singleHap[i].orUp) + "]";
    data.push_back(OR);

    if (this->adjust) {
      data.push_back(convert2string(this->Results.singleHap[i].HolmP));
      data.push_back(convert2string(this->Results.singleHap[i].SidakSSP));
      data.push_back(convert2string(this->Results.singleHap[i].SidakSDP));
      data.push_back(convert2string(this->Results.singleHap[i].BHP));
      data.push_back(convert2string(this->Results.singleHap[i].BYP));
    }

    html->addDataRow(data);
  }
  return html->getTable();
}

std::string HaplotypeBase::reporttxttableQTL() {
  std::stringstream res;
  res.precision(3);
  res << "Haplotype\tTotal count\tRegression coefficient\tStandard error\t "
         "Regression r2\t T statistics\tP value";
  if (this->adjust) {
    res << "\tHolm\tSidakSS\tSidakSD\tFDR_BH\tFDR_BY";
  }
  res << "\n";
  for (int i = 0; i < this->Results.singleHapQTL.size(); i++) {
    if (0 == this->Results.singleHapQTL[i].ValidSampleNum) continue;
    for (int j = 0; j < this->SnpIdx.size(); j++) {
      res << this->data->getallele(this->Results.haplotypes[i][j]);
    }
    res << "\t\t";
    res << this->Results.BothCount[i] << "\t\t";
    res << this->Results.singleHapQTL[i].beta << "\t\t\t";
    res << this->Results.singleHapQTL[i].se << "\t\t";
    res << this->Results.singleHapQTL[i].R2 << "\t\t";
    res << this->Results.singleHapQTL[i].T << "\t\t";
    res << this->Results.singleHapQTL[i].p;
    if (this->adjust) {
      res << "\t" << this->Results.singleHapQTL[i].HolmP << "\t"
          << this->Results.singleHapQTL[i].SidakSSP << "\t"
          << this->Results.singleHapQTL[i].SidakSDP << "\t"
          << this->Results.singleHapQTL[i].BHP << "\t"
          << this->Results.singleHapQTL[i].BYP;
    }
    res << "\n";
  };
  return res.str();
}

std::string HaplotypeBase::reporthtmltableQTL() {
  boost::shared_ptr<SHEsis::CreatHtmlTable> html(new SHEsis::CreatHtmlTable());
  html->createTable("Haplotype_Analysis");
  std::vector<std::string> data;
  data.push_back("Haplotype");
  data.push_back("Total count");
  data.push_back("Beta");
  data.push_back("SE");
  data.push_back("R<sup>2</sup>");
  data.push_back("T");
  data.push_back("p");
  if (this->adjust) {
    data.push_back("Holm");
    data.push_back("SidakSS");
    data.push_back("SidakSD");
    data.push_back("FDR_BH");
    data.push_back("FDR_BY");
  }
  html->addHeadRow(data);
  for (int i = 0; i < this->Results.singleHapQTL.size(); i++) {
    if (0 == this->Results.singleHapQTL[i].ValidSampleNum) continue;
    data.clear();
    std::stringstream hap;
    for (int j = 0; j < this->SnpIdx.size(); j++) {
      hap << this->data->getallele(this->Results.haplotypes[i][j]);
    }
    data.push_back(hap.str());
    data.push_back(convert2string(this->Results.BothCount[i]));
    data.push_back(convert2string(this->Results.singleHapQTL[i].beta));
    data.push_back(convert2string(this->Results.singleHapQTL[i].se));
    data.push_back(convert2string(this->Results.singleHapQTL[i].R2));
    data.push_back(convert2string(this->Results.singleHapQTL[i].T));
    data.push_back(convert2string(this->Results.singleHapQTL[i].p));
    if (this->adjust) {
      data.push_back(convert2string(this->Results.singleHapQTL[i].HolmP));
      data.push_back(convert2string(this->Results.singleHapQTL[i].SidakSSP));
      data.push_back(convert2string(this->Results.singleHapQTL[i].SidakSDP));
      data.push_back(convert2string(this->Results.singleHapQTL[i].BHP));
      data.push_back(convert2string(this->Results.singleHapQTL[i].BYP));
    }
    html->addDataRow(data);
  }
  return html->getTable();
}
}
