/*
 * SHEsisData.cpp
 *
 *  Created on: Aug 3, 2014
 *      Author: ionadmin
 */

#include "SHEsisData.h"
#include <algorithm>
#include <boost/assert.hpp>
#include <sstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
namespace SHEsis {

SHEsisData::SHEsisData(int SampleNum, int SnpNum, int NumOfChrSet)
    : SampleNum(SampleNum),
      SnpNum(SnpNum),
      NumOfChrSet(NumOfChrSet),
      mGenotype(boost::extents[SampleNum][SnpNum][NumOfChrSet]),
      vLocusInfo(SnpNum),
      CaseNum(-1),
      ControlNum(-1),
      codeIdx(1),
      vLabel(SampleNum),
      vLocusName(SnpNum),
      missingalleles(SnpNum) {
  for (int i = 0; i < this->vLocusInfo.size(); i++) {
    std::stringstream ss;
    ss << "site" << i + 1;
    this->vLocusName[i] = ss.str();
  }
};

SHEsisData::~SHEsisData() {
  this->vLabel.clear();
  this->vLocusInfo.clear();
}
void SHEsisData::setLocusName(int snpidx, std::string s) {
  this->vLocusName[snpidx] = s;
}
std::string SHEsisData::getOriginGenotype(std::string geno) {
  std::vector<std::string> strs;
  std::string res;
  boost::split(strs, geno, boost::is_any_of("/"));
  for (int i = 0; i < strs.size() - 1; i++)
    res += this->getallele(std::atoi(strs[i].c_str())) + "/";
  res += this->getallele(std::atoi(strs[strs.size() - 1].c_str()));
  return res;
}

std::string SHEsisData::getStrFromSortedGenotype(std::vector<short> v) {
  std::stringstream ss;
  for (int i = 0; i < v.size() - 1; i++) {
    ss << v[i] << "/";
  }
  ss << v[v.size() - 1];
  return ss.str();
}

short SHEsisData::GetAlleleCode(std::string const val) {
  if (std::strcmp("0", val.c_str()) == 0) return 0;
  boost::unordered_map<short, std::string>::iterator iter;
  for (iter = this->code2allele.begin(); iter != this->code2allele.end();
       iter++) {
    if (val == iter->second) return iter->first;
  }
  this->code2allele[this->codeIdx] = val;
  return this->codeIdx++;
}

void SHEsisData::printLocusInfo() {
  for (int i = 0; i < this->vLocusInfo.size(); i++) {
    boost::unordered_map<short, double>::iterator map_it;
    std::cout << "\nLocus " << i;
    std::cout << "\nCaseAlleleCount:\n";
    for (map_it = this->vLocusInfo[i].CaseAlleleCount.begin();
         map_it != this->vLocusInfo[i].CaseAlleleCount.end(); map_it++) {
      std::cout << map_it->first << ":" << map_it->second << ", ";
    }
    std::cout << "\nControlAlleleCount:\n";
    for (map_it = this->vLocusInfo[i].ControlAlleleCount.begin();
         map_it != this->vLocusInfo[i].ControlAlleleCount.end(); map_it++) {
      std::cout << map_it->first << ":" << map_it->second << ", ";
    }
    std::cout << "\nBothAlleleCount:\n";
    for (map_it = this->vLocusInfo[i].BothAlleleCount.begin();
         map_it != this->vLocusInfo[i].BothAlleleCount.end(); map_it++) {
      std::cout << map_it->first << ":" << map_it->second << ", ";
    }
    boost::unordered_map<std::string, double>::iterator map_it2;
    std::cout << "\nCaseGenotypeCount:\n";
    for (map_it2 = this->vLocusInfo[i].CaseGenotypeCount.begin();
         map_it2 != this->vLocusInfo[i].CaseGenotypeCount.end(); map_it2++) {
      std::cout << map_it2->first << ":" << map_it2->second << ", ";
    }
    std::cout << "\nControlGenotypeCount:\n";
    for (map_it2 = this->vLocusInfo[i].ControlGenotypeCount.begin();
         map_it2 != this->vLocusInfo[i].ControlGenotypeCount.end(); map_it2++) {
      std::cout << map_it2->first << ":" << map_it2->second << ", ";
    }
    std::cout << "\nBothGenotypeCount:\n";
    for (map_it2 = this->vLocusInfo[i].BothGenotypeCount.begin();
         map_it2 != this->vLocusInfo[i].BothGenotypeCount.end(); map_it2++) {
      std::cout << map_it2->first << ":" << map_it2->second << ", ";
    }
  }
}

void SHEsisData::getCaseAndControlNum() {
  this->CaseNum = 0;
  this->ControlNum = 0;
  for (int iSample = 0; iSample < this->SampleNum; iSample++) {
    if (CASE == this->vLabel[iSample]) {
      this->CaseNum++;
    } else if (CONTROL == this->vLabel[iSample]) {
      this->ControlNum++;
    }
  }
}

int SHEsisData::getCaseNum() {
  if (-1 == this->CaseNum) this->getCaseAndControlNum();
  return this->CaseNum;
}

double SHEsisData::getCallrate(int snp) {
  return 1 - (double)(this->missingalleles[snp].CaseAlleleNum +
                      this->missingalleles[snp].CtrlAlleleNum) /
                 (this->getSampleNum() * this->getNumOfChrSet());
  //	return
  //1-(double)(this->missingalleles[snp].CaseAlleleNum)/(double)(this->getSampleNum()*this->getNumOfChrSet());
}

double SHEsisData::getCaseCallrate(int snp) {
  return 1 - (double)(this->missingalleles[snp].CaseAlleleNum) /
                 (double)(this->getCaseNum() * this->getNumOfChrSet());
}

int SHEsisData::getControlNum() {
  if (-1 == this->ControlNum) this->getCaseAndControlNum();
  return this->ControlNum;
}

double SHEsisData::getControlCallrate(int snp) {
  return 1 - (double)(this->missingalleles[snp].CtrlAlleleNum) /
                 (double)(this->getControlNum() * this->getNumOfChrSet());
}

void SHEsisData::statCount() {
  for (int iSnp = 0; iSnp < this->SnpNum; iSnp++) {
    this->missingalleles[iSnp].CaseAlleleNum = 0;
    this->missingalleles[iSnp].CtrlAlleleNum = 0;
    this->vLocusInfo[iSnp].AlleleCallrate = 0;
  }
  this->vLocusInfo.clear();
  this->vLocusInfo.resize(SnpNum);
  std::vector<short> geno;
  for (int iSample = 0; iSample < this->SampleNum; iSample++) {
    for (int iSnp = 0; iSnp < this->SnpNum; iSnp++) {
      geno.clear();
      for (int iChrset = 0; iChrset < this->NumOfChrSet; iChrset++) {
        BOOST_ASSERT(iSample < this->mGenotype.shape()[0] &&
                     iSnp < this->mGenotype.shape()[1] &&
                     iChrset < this->mGenotype.shape()[2]);
        short Cur = this->mGenotype[iSample][iSnp][iChrset];
        if (GENOTYPE_MISSING == Cur) {
          this->vLocusInfo[iSnp].AlleleCallrate++;
          this->missingalleles[iSnp].CaseAlleleNum++;  // for qtl, using
                                                       // CaseAlleleNum to store
                                                       // the missingness status
          continue;
        }
        {
          if (this->vLocusInfo[iSnp].BothAlleleCount.end() ==
              this->vLocusInfo[iSnp].BothAlleleCount.find(Cur))
            this->vLocusInfo[iSnp].BothAlleleCount[Cur] = 1;
          else
            (this->vLocusInfo[iSnp].BothAlleleCount[Cur])++;

          if (iChrset == (this->NumOfChrSet - 1)) {
            geno.push_back(Cur);
            if (this->NumOfChrSet > geno.size()) {
              this->vLocusInfo[iSnp].GenoCallrate++;
              continue;  // indicating there is missing data within this site
                         // for the specific individual
            }
            std::sort(geno.begin(), geno.end());
            std::string genoStr = getStrFromSortedGenotype(geno);
            if (this->vLocusInfo[iSnp].BothGenotypeCount.end() ==
                this->vLocusInfo[iSnp].BothGenotypeCount.find(genoStr))
              this->vLocusInfo[iSnp].BothGenotypeCount[genoStr] = 1;
            else
              (this->vLocusInfo[iSnp].BothGenotypeCount[genoStr])++;
          } else
            geno.push_back(Cur);
        }
      }
    }
  }
  for (int iSnp = 0; iSnp < this->SnpNum; iSnp++) {
    this->vLocusInfo[iSnp].AlleleCallrate =
        1 - (this->vLocusInfo[iSnp].AlleleCallrate /
             (double)this->getNumOfChrSet() / this->getSampleNum());
    this->vLocusInfo[iSnp].GenoCallrate =
        1 -
        (this->vLocusInfo[iSnp].GenoCallrate / (double)this->getSampleNum());
  }
}

void SHEsisData::statCount(std::vector<SampleStatus>& label) {
  for (int iSnp = 0; iSnp < this->SnpNum; iSnp++) {
    this->missingalleles[iSnp].CaseAlleleNum = 0;
    this->missingalleles[iSnp].CtrlAlleleNum = 0;
    this->vLocusInfo[iSnp].AlleleCallrate = 0;
  }
  this->vLocusInfo.clear();
  this->vLocusInfo.resize(SnpNum);
  std::vector<short> geno;
  for (int iSample = 0; iSample < this->SampleNum; iSample++) {
    for (int iSnp = 0; iSnp < this->SnpNum; iSnp++) {
      geno.clear();
      for (int iChrset = 0; iChrset < this->NumOfChrSet; iChrset++) {
        BOOST_ASSERT(iSample < this->mGenotype.shape()[0] &&
                     iSnp < this->mGenotype.shape()[1] &&
                     iChrset < this->mGenotype.shape()[2]);
        short Cur = this->mGenotype[iSample][iSnp][iChrset];
        if (GENOTYPE_MISSING == Cur) {
          this->vLocusInfo[iSnp].AlleleCallrate++;
          if (CASE == label[iSample])
            this->missingalleles[iSnp].CaseAlleleNum++;
          else if (CONTROL == label[iSample])
            this->missingalleles[iSnp].CtrlAlleleNum++;
          continue;
        }
        if (CASE == label[iSample]) {
          if (this->vLocusInfo[iSnp].CaseAlleleCount.end() ==
              this->vLocusInfo[iSnp].CaseAlleleCount.find(Cur))
            this->vLocusInfo[iSnp].CaseAlleleCount[Cur] = 1;
          else
            (this->vLocusInfo[iSnp].CaseAlleleCount[Cur])++;
          if (this->vLocusInfo[iSnp].ControlAlleleCount.end() ==
              this->vLocusInfo[iSnp].ControlAlleleCount.find(Cur))
            this->vLocusInfo[iSnp].ControlAlleleCount[Cur] = 0;
          if (iChrset == (this->NumOfChrSet - 1)) {
            geno.push_back(Cur);
            if (this->NumOfChrSet > geno.size()) {
              this->vLocusInfo[iSnp].GenoCallrate++;
              continue;  // indicating there is missing data within this site
                         // for the specific individual
            }

            std::sort(geno.begin(), geno.end());
            std::string genoStr = getStrFromSortedGenotype(geno);
            if (this->vLocusInfo[iSnp].CaseGenotypeCount.end() ==
                this->vLocusInfo[iSnp].CaseGenotypeCount.find(genoStr))
              this->vLocusInfo[iSnp].CaseGenotypeCount[genoStr] = 1;
            else
              (this->vLocusInfo[iSnp].CaseGenotypeCount[genoStr])++;

            if (this->vLocusInfo[iSnp].ControlGenotypeCount.end() ==
                this->vLocusInfo[iSnp].ControlGenotypeCount.find(genoStr))
              this->vLocusInfo[iSnp].ControlGenotypeCount[genoStr] = 0;
          } else
            geno.push_back(Cur);

          continue;
        };
        if (CONTROL == label[iSample]) {
          if (this->vLocusInfo[iSnp].ControlAlleleCount.end() ==
              this->vLocusInfo[iSnp].ControlAlleleCount.find(Cur))
            this->vLocusInfo[iSnp].ControlAlleleCount[Cur] = 1;
          else
            (this->vLocusInfo[iSnp].ControlAlleleCount[Cur])++;
          if (this->vLocusInfo[iSnp].CaseAlleleCount.end() ==
              this->vLocusInfo[iSnp].CaseAlleleCount.find(Cur))
            this->vLocusInfo[iSnp].CaseAlleleCount[Cur] = 0;
          if (iChrset == (this->NumOfChrSet - 1)) {
            geno.push_back(Cur);
            if (this->NumOfChrSet > geno.size()) {
              this->vLocusInfo[iSnp].GenoCallrate++;
              continue;  // indicating there is missing data within this site
                         // for the specific individual
            }

            std::sort(geno.begin(), geno.end());
            std::string genoStr = getStrFromSortedGenotype(geno);
            if (this->vLocusInfo[iSnp].ControlGenotypeCount.end() ==
                this->vLocusInfo[iSnp].ControlGenotypeCount.find(genoStr))
              this->vLocusInfo[iSnp].ControlGenotypeCount[genoStr] = 1;
            else
              (this->vLocusInfo[iSnp].ControlGenotypeCount[genoStr])++;
            if (this->vLocusInfo[iSnp].CaseGenotypeCount.end() ==
                this->vLocusInfo[iSnp].CaseGenotypeCount.find(genoStr))
              this->vLocusInfo[iSnp].CaseGenotypeCount[genoStr] = 0;
          } else
            geno.push_back(Cur);
        }
      }
    }
  }
  for (int iSnp = 0; iSnp < this->SnpNum; iSnp++) {
    BOOST_ASSERT(this->vLocusInfo[iSnp].CaseGenotypeCount.size() ==
                 this->vLocusInfo[iSnp].ControlGenotypeCount.size());
    boost::unordered_map<std::string, double>::iterator iter;
    for (iter = this->vLocusInfo[iSnp].CaseGenotypeCount.begin();
         iter != this->vLocusInfo[iSnp].CaseGenotypeCount.end(); iter++) {
      this->vLocusInfo[iSnp].BothGenotypeCount[iter->first] =
          this->vLocusInfo[iSnp].CaseGenotypeCount[iter->first] +
          this->vLocusInfo[iSnp].ControlGenotypeCount[iter->first];
    }

    BOOST_ASSERT(this->vLocusInfo[iSnp].CaseAlleleCount.size() ==
                 this->vLocusInfo[iSnp].ControlAlleleCount.size());
    boost::unordered_map<short, double>::iterator iter2;
    for (iter2 = this->vLocusInfo[iSnp].CaseAlleleCount.begin();
         iter2 != this->vLocusInfo[iSnp].CaseAlleleCount.end(); iter2++) {
      this->vLocusInfo[iSnp].BothAlleleCount[iter2->first] =
          this->vLocusInfo[iSnp].CaseAlleleCount[iter2->first] +
          this->vLocusInfo[iSnp].ControlAlleleCount[iter2->first];
    }
  }

  for (int iSnp = 0; iSnp < this->SnpNum; iSnp++) {
    this->vLocusInfo[iSnp].AlleleCallrate =
        1 - (this->vLocusInfo[iSnp].AlleleCallrate /
             (double)this->getNumOfChrSet() / this->getSampleNum());
    this->vLocusInfo[iSnp].GenoCallrate =
        1 -
        (this->vLocusInfo[iSnp].GenoCallrate / (double)this->getSampleNum());
  }
}
}
