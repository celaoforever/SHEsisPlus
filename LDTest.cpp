/*
 * LDTest.cpp
 *
 *  Created on: Sep 17, 2014
 *      Author: ada
 */

#include "LDTest.h"
//#include "HaplotypeDiploid.h"
#include "HaplotypeEM.h"
#include "Haplotype.h"
#include <boost/shared_ptr.hpp>
#include "CreatHtmlTable.h"
#define PI 3.1415926535898
#define ABS(x) x < 0 ? (-1 * (x)) : (x)
namespace SHEsis {

LDTest::LDTest(boost::shared_ptr<SHEsisData> data, std::string path)
    : data(data),
      ldtype(LD_IN_BOTH),
      path(path),
      bForceSAT(false),
      resD(boost::extents[this->data->getSnpNum()][this->data->getSnpNum()]),
      resR2(boost::extents[this->data->getSnpNum()][this->data->getSnpNum()]) {
  // TODO Auto-generated constructor stub
}

LDTest::~LDTest() {
  // TODO Auto-generated destructor stub
}

void LDTest::printRes() {

  std::cout << "D':\n";
  for (int i = 0; i < this->resD.shape()[0]; i++) {
    for (int j = 0; j < this->resD.shape()[1]; j++) {
      std::cout << this->resD[i][j] << "\t";
    }
    std::cout << "\n";
  }
  std::cout << "R2:\n";
  for (int i = 0; i < this->resR2.shape()[0]; i++) {
    for (int j = 0; j < this->resR2.shape()[1]; j++) {
      std::cout << this->resR2[i][j] << "\t";
    }
    std::cout << "\n";
  }
}

std::string LDTest::reporttxt() {
  std::stringstream res;
  res.precision(3);
  res << "\n-------------------------------------------------------\n";
  res << "Linkage Disequilibrium Analysis\n";
  res << "-------------------------------------------------------\n";
  res << "D'\t";
  for (int i = 0; i < this->data->vLocusName.size(); i++) {
    res << this->data->vLocusName[i] << "\t";
  }
  res << "\n";
  for (int i = 0; i < this->resD.shape()[0]; i++) {
    res << this->data->vLocusName[i] << "\t";
    for (int j = 0; j < this->resD.shape()[1]; j++) {
      if (i < j)
        res << this->resD[i][j] << "\t";
      else
        res << "-"
            << "\t";
    }
    res << "\n";
  }
  res << "-------------------------------------------------------\n";
  res << "R2\t";
  for (int i = 0; i < this->data->vLocusName.size(); i++) {
    res << this->data->vLocusName[i] << "\t";
  }
  res << "\n";
  for (int i = 0; i < this->resR2.shape()[0]; i++) {
    res << this->data->vLocusName[i] << "\t";
    for (int j = 0; j < this->resR2.shape()[1]; j++) {
      if (i < j)
        res << this->resR2[i][j] << "\t";
      else
        res << "-"
            << "\t";
    }
    res << "\n";
  }
  res << "-------------------------------------------------------\n";
  return res.str();
}

int getHapIdx(std::vector<boost::shared_ptr<short[]> >& v, int a1, int a2) {
  for (int i = 0; i < v.size(); i++) {
    if (a1 == v[i][0] && a2 == v[i][1]) return i;
  }
  return -1;
}

void LDTest::resetMissing() {
  this->Snp1AlleleCount.casecount.clear();
  this->Snp1AlleleCount.ctrlcount.clear();
  this->Snp2AlleleCount.casecount.clear();
  this->Snp2AlleleCount.ctrlcount.clear();
};

void LDTest::StatAllelesInNonmissingSample(boost::shared_ptr<bool[]> missing,
                                           int snp1, int snp2) {
  this->resetMissing();
  int sampleNum = this->data->getSampleNum();
  int ploidy = this->data->getNumOfChrSet();
  for (int iSample = 0; iSample < sampleNum; iSample++) {
    if (missing[iSample]) continue;
    for (int p = 0; p < ploidy; p++) {
      short curAllele1 = this->data->mGenotype[iSample][snp1][p];
      short curAllele2 = this->data->mGenotype[iSample][snp2][p];
      if (this->data->vQuantitativeTrait.size() > 0 ||
          CASE == this->data->vLabel[iSample]) {
        if (this->Snp1AlleleCount.casecount.find(curAllele1) ==
            this->Snp1AlleleCount.casecount.end())
          this->Snp1AlleleCount.casecount[curAllele1] = 1;
        else
          this->Snp1AlleleCount.casecount[curAllele1]++;
        if (this->Snp2AlleleCount.casecount.find(curAllele2) ==
            this->Snp2AlleleCount.casecount.end())
          this->Snp2AlleleCount.casecount[curAllele2] = 1;
        else
          this->Snp2AlleleCount.casecount[curAllele2]++;
      } else if (CONTROL == this->data->vLabel[iSample]) {
        if (this->Snp1AlleleCount.ctrlcount.find(curAllele1) ==
            this->Snp1AlleleCount.ctrlcount.end())
          this->Snp1AlleleCount.ctrlcount[curAllele1] = 1;
        else
          this->Snp1AlleleCount.ctrlcount[curAllele1]++;
        if (this->Snp2AlleleCount.ctrlcount.find(curAllele2) ==
            this->Snp2AlleleCount.ctrlcount.end())
          this->Snp2AlleleCount.ctrlcount[curAllele2] = 1;
        else
          this->Snp2AlleleCount.ctrlcount[curAllele2]++;
      }
    }
  }
}

void LDTest::TwoLociLDTest(int snp1, int snp2, LD_TYPE type, double& R2,
                           double& D) {
  //	std::cout<<"site "<<snp1<<"-site"<<snp2<<"\n";
  D = 0;
  R2 = 0;
  std::vector<short> mask;
  for (int i = 0; i < this->data->getSnpNum(); i++) {
    if (i == snp1 || i == snp2)
      mask.push_back(1);
    else
      mask.push_back(0);
  }
  BOOST_ASSERT(mask.size() == this->data->getSnpNum());
  if (this->bForceSAT) {
    this->hp.reset(new Haplotype(this->data, 2, mask));
  } else {
    this->hp.reset(new HaplotypeEM(this->data, 2, mask));
  }
  //  if (this->data->getNumOfChrSet() <= 2) {
  //    if (this->bForceSAT)
  //      this->hp.reset(new Haplotype(this->data, 2, mask));
  //    else
  //      this->hp.reset(new HaplotypeDiploid(this->data, 2, mask));
  //  } else {
  //    this->hp.reset(new Haplotype(this->data, 2, mask));
  //  }
  hp->startHaplotypeAnalysis();
  boost::shared_ptr<HaplotypeEM> hp_em =
      boost::dynamic_pointer_cast<HaplotypeEM>(this->hp);
  if (hp_em) StatAllelesInNonmissingSample(hp_em->missing, snp1, snp2);

  boost::unordered_map<short, double>::iterator iter1;
  boost::unordered_map<short, double>::iterator iter2;

  for (iter1 = this->data->vLocusInfo[snp1].BothAlleleCount.begin();
       iter1 != this->data->vLocusInfo[snp1].BothAlleleCount.end(); iter1++) {
    for (iter2 = this->data->vLocusInfo[snp2].BothAlleleCount.begin();
         iter2 != this->data->vLocusInfo[snp2].BothAlleleCount.end(); iter2++) {
      int i = getHapIdx(hp->Results.haplotypes, iter1->first, iter2->first);
      double hapfreq = 0;
      double allelefreq1 = 0;
      double allelefreq2 = 0;
      int allele = 0;
      double d = 0;
      double dmax = 0;
      switch (type) {
        case LD_IN_CASE:

          hapfreq = i == -1 ? 0 : (double)hp->Results.CaseCount[i] /
                                      (double)(this->data->getNumOfChrSet() *
                                               /*this->data->getCaseNum()*/ this
                                                   ->hp->NonmissingCase);
          allele = iter1->first;
          allelefreq1 =
              hp_em ? ((double)this->Snp1AlleleCount.casecount[allele] /
                       (double)(this->hp->NonmissingCase *
                                this->data->getNumOfChrSet()))
                    : ((double)data->vLocusInfo[snp1].CaseAlleleCount[allele] /
                       (double)(this->data->getNumOfChrSet() *
                                this->data->getCaseCallrate(snp1) *
                                this->data->getCaseNum()));
          allele = iter2->first;
          allelefreq2 =
              hp_em ? ((double)this->Snp2AlleleCount.casecount[allele] /
                       (double)(this->hp->NonmissingCase *
                                this->data->getNumOfChrSet()))
                    : ((double)data->vLocusInfo[snp2].CaseAlleleCount[allele] /
                       ((double)(this->data->getNumOfChrSet() *
                                 this->data->getCaseCallrate(snp2) *
                                 this->data->getCaseNum())));
          break;

        case LD_IN_CTRL:
          hapfreq =
              i == -1 ? 0 : (double)hp->Results.ControlCount[i] /
                                (double)(this->data->getNumOfChrSet() *
                                         /*this->data->getControlNum()*/ this
                                             ->hp->NonmissingCtrl);
          allele = iter1->first;
          allelefreq1 =
              hp_em
                  ? ((double)this->Snp1AlleleCount.ctrlcount[allele] /
                     (double)(this->hp->NonmissingCtrl *
                              this->data->getNumOfChrSet()))
                  : ((double)data->vLocusInfo[snp1].ControlAlleleCount[allele] /
                     (double)(this->data->getNumOfChrSet() *
                              this->data->getControlCallrate(snp1) *
                              this->data->getControlNum()));
          allele = iter2->first;
          allelefreq2 =
              hp_em
                  ? ((double)this->Snp2AlleleCount.ctrlcount[allele] /
                     (double)(this->hp->NonmissingCtrl *
                              this->data->getNumOfChrSet()))
                  : ((double)data->vLocusInfo[snp2].ControlAlleleCount[allele] /
                     (double)(this->data->getNumOfChrSet() *
                              this->data->getControlCallrate(snp2) *
                              this->data->getControlNum()));
          break;
        case LD_IN_BOTH:
          hapfreq = i == -1 ? 0 : (double)(hp->Results.BothCount[i]) /
                                      (double)(this->data->getNumOfChrSet() *
                                               /*this->data->getSampleNum()*/
                                               (this->hp->NonmissingCase +
                                                this->hp->NonmissingCtrl));
          allele = iter1->first;
          allelefreq1 =
              hp_em ? (this->Snp1AlleleCount.casecount[allele] +
                       this->Snp1AlleleCount.ctrlcount[allele]) /
                          (double)((this->hp->NonmissingCtrl +
                                    this->hp->NonmissingCase) *
                                   this->data->getNumOfChrSet())
                    : ((double)data->vLocusInfo[snp1].BothAlleleCount[allele] /
                       (double)(this->data->getNumOfChrSet() *
                                this->data->getCallrate(snp1) *
                                this->data->getSampleNum()));
          allele = iter2->first;
          allelefreq2 =
              hp_em ? (this->Snp2AlleleCount.casecount[allele] +
                       this->Snp2AlleleCount.ctrlcount[allele]) /
                          (double)((this->hp->NonmissingCtrl +
                                    this->hp->NonmissingCase) *
                                   this->data->getNumOfChrSet())
                    : (double)data->vLocusInfo[snp2].BothAlleleCount[allele] /
                          (double)(this->data->getNumOfChrSet() *
                                   this->data->getCallrate(snp2) *
                                   this->data->getSampleNum());
          //          std::cout<<"hapfreq="<<hapfreq<<",a1="<<allelefreq1<<",a2="<<allelefreq2;
          break;
        default:
          std::cout << "***ERROR:no such ld analysis type\n";
          exit(-1);
      };
      //      std::cout<<"Nonmissing
      // case="<<this->hp->NonmissingCase<<",Nonmissing
      // ctrl="<<this->hp->NonmissingCtrl<<",allelfreq1="<<allelefreq1<<","<<"allelefreq2="<<allelefreq2<<",hapfreq="<<hapfreq<<"\n";
      d = hapfreq - allelefreq1 * allelefreq2;
      double pq = allelefreq1 * allelefreq2;
      double p1q1 = (1 - allelefreq1) * (1 - allelefreq2);
      double p1q = allelefreq1 * (1 - allelefreq2);
      double q1p = allelefreq2 * (1 - allelefreq1);
      if (d < 0) {
        dmax = pq < p1q1 ? pq : p1q1;
      } else if (d > 0) {
        dmax = p1q < q1p ? p1q : q1p;
      }
      if (dmax != 0) {
        double absd = ABS(d / dmax);
        D += allelefreq1 * allelefreq2 * absd;
        R2 += d * d / (1 - allelefreq1) / (1 - allelefreq2);
      }
      //      std::cout<<",d="<<d<<",dmax="<<dmax<<",d/dmax="<<adsd<<",normalizedd="<<normalizedD<<"\n";
    }
  }
}

void LDTest::AllLociLDtest() {
  if (this->data->vLocusInfo[0].BothAlleleCount.size() == 0 &&
      this->data->vQuantitativeTrait.size() > 0) {
    this->data->statCount();
    this->ldtype = LD_IN_BOTH;
  } else if (this->data->vLocusInfo[0].BothAlleleCount.size() == 0 &&
             this->data->vQuantitativeTrait.size() == 0)
    this->data->statCount(this->data->vLabel);

  //  if(data->vQuantitativeTrait.size()>0){
  //
  //	  this->data->statCount();
  //	  this->ldtype=LD_IN_BOTH;
  //  }
  //  else
  //	  this->data->statCount(this->data->vLabel);
  int total = this->data->getSnpNum() * (this->data->getSnpNum() - 1) / 2;
  int count = 0;
  for (int i = 0; i < this->data->getSnpNum(); i++) {
    for (int j = i + 1; j < this->data->getSnpNum(); j++) {
      double D, R2;
      this->TwoLociLDTest(i, j, this->ldtype, R2, D);
      this->resD[i][j] = D;
      this->resR2[i][j] = R2;
      printf("\rProgress:%d%%", (int)(100 * (double)count / (double)total));
      fflush(stdout);
      count++;
    }
  }
  printf("\rProgress:%d%%\n", 100);
  fflush(stdout);
}
bool existsNonZeroPixel(int x_start, int x_end, int y, int length,
                        boost::shared_ptr<bool[]> p) {
  for (int x = x_start; x <= x_end; x++) {
    if (p[y * length + x]) return true;
  }
  return false;
}
void DrawSquare(BMP* bmp, RGB rgb, int center_x, int center_y, double l) {
  double length = ceil(sqrt(2) * l + 5);
  boost::shared_ptr<bool[]> p(new bool[(int)length * (int)length]);
  boost::shared_ptr<bool[]> p_rotated(new bool[(int)length * (int)length]);
  for (int i = 0; i < length * length; i++) {
    p[i] = false;
    p_rotated[i] = false;
  }
  double cx = length / 2;
  double cy = cx;

  // draw square
  int x_start = (int)round(cx - (double)l / 2.0);
  int x_end = (int)round(cx + (double)l / 2.0);
  int y = x_start;
  for (int x = x_start; x < x_end; x++) {
    p[length * y + x] = true;
  }
  y = y + l;
  for (int x = x_start; x < x_end; x++) {
    p[length * y + x] = true;
  }
  int y_start = x_start;
  int y_end = y_start + l;
  int x = y_start;
  for (int y = y_start; y < y_end; y++) {
    p[length * y + x] = true;
  };
  x = x + l;
  for (int y = y_start; y < y_end; y++) {
    p[length * y + x] = true;
  };

  // rotate
  double r = 45.0 / 180.0 * PI;
  int ymax = 0;
  int ymax_x = 0;
  int ymin = length;
  int ymin_x = 0;
  for (int x = 0; x < length; x++) {
    for (int y = 0; y < length; y++)
      if (p[y * length + x]) {
        int newx = (int)round((x - cx) * cos(r) - (y - cy) * sin(r) + cx);
        int newy = (int)round((x - cx) * sin(r) + (y - cy) * cos(r) + cy);
        if (ymax < newy) {
          ymax = newy;
          ymax_x = newx;
        };
        if (ymin > newy) {
          ymin = newy;
          ymin_x = newx;
        }
        p_rotated[newy * length + newx] = true;
      }
  }

  // fill
  for (int y = ymin + 1; y < ymax; y++) {
    for (x = ymin_x; x >= 0; x--) {
      if (y < ymin + 3 &&
          !existsNonZeroPixel(x - 10, x, y, length, p_rotated)) {
        break;
      }
      if (y > ymax - 3 &&
          !existsNonZeroPixel(x - 10, x, y, length, p_rotated)) {
        break;
      }
      if (!p_rotated[y * length + x])
        p_rotated[y * length + x] = true;
      else
        break;
    }

    for (x = ymin_x + 1; x < length; x++) {
      if (y > ymax - 3 &&
          !existsNonZeroPixel(x, x + 10, y, length, p_rotated)) {
        break;
      }
      if (y < ymin + 3 &&
          !existsNonZeroPixel(x, x + 10, y, length, p_rotated)) {
        break;
      }
      if (!p_rotated[y * length + x])
        p_rotated[y * length + x] = true;
      else
        break;
    }
  };

  // draw on bmp
  int translate_x = center_x - cx;
  int translate_y = center_y - cy;
  for (int x = 0; x < length; x++) {
    for (int y = 0; y < length; y++) {
      if (p_rotated[length * y + x]) {
        BMP_point(bmp, x + translate_x, y + translate_y, rgb);
      }
    }
  }
}

void LDTest::DrawLDMap(int type) {  // 0 for D', 1 for R2
  double snpnum = (double)this->data->getSnpNum();
  double height, width;
  int recnum;
  recnum = snpnum - 1;
  if (snpnum < 10)
    height = 600;
  else
    height = snpnum * 60;
  width = (height - 230) / snpnum * 2 * (snpnum - 1) + 100;
  double sidelength = (height - 230) / snpnum * 2;
  double insidelength = sidelength - 2;
  double from_x = (double)width / 2 - (double)recnum * (double)sidelength / 2;
  this->ldmap = BMP_new(width, height);
  BMP_clear(this->ldmap, 0xeeeeee);

  for (int i = 0; i < snpnum; i++) {
    int step = i * sidelength;
    std::stringstream ss;
    ss << (i + 1);
    BMP_draw_string(this->ldmap, ss.str().c_str(), from_x + step - 4, 200 - 11,
                    RGB_BLACK, 0);
    BMP_line(this->ldmap, from_x + step, 20, from_x + step, 25, RGB_BLACK);
    int width = 115 - this->data->vLocusName[i].length() * 5;
    width = width > 0 ? width : 0;
    BMP_line(this->ldmap, from_x + step, 25, from_x + step, width, RGB_GRAY);
    BMP_draw_string(this->ldmap, this->data->vLocusName[i].c_str(),
                    from_x + step - 7, 140, RGB_BLACK, -90.0 / 180.0 * PI);
  }
  BMP_line(this->ldmap, from_x, 20, from_x + recnum * sidelength, 20, RGB_GRAY);
  BMP_line(this->ldmap, from_x, 25, from_x + recnum * sidelength, 25, RGB_GRAY);

  for (int j = 0; j < recnum; j++) {
    int count = 0;
    for (int i = 1; i <= recnum - j; i++) {
      int step = i * sidelength;
      int begin_x = from_x + sidelength / 2 * j;
      int begin_y = 200 + sidelength / 2 * j;
      double score = type == 0 ? this->resD[count][count + j + 1]
                               : this->resR2[count][count + j + 1];
      RGB color = (1 - score) * 0xff;
      color = color << 8;
      color += (1 - score) * 0xff + 0xff0000;
      DrawSquare(this->ldmap, color, begin_x + step - sidelength / 2,
                 begin_y + sidelength / 2, sidelength / sqrt(2) + 2);
      std::stringstream ss;
      double s = (int)(score * 100);
      ss << (s / 100.0);
      int strx = begin_x + step - sidelength / 2 - 7;
      int stry = begin_y + sidelength / 2 - 7;
      BMP_draw_string(this->ldmap, ss.str().c_str(), strx, stry, RGB_BLACK, 0);
      count++;
    }
  }
  if (0 == type) {
    BMP_draw_string(this->ldmap, "D'", 30, height - 50, RGB_BLACK, 0);
  } else if (1 == type) {
    BMP_draw_string(this->ldmap, "R", 30, height - 50, RGB_BLACK, 0);
    BMP_draw_mini_string(this->ldmap, "2", 45, height - 55, RGB_BLACK);
  } else
    BOOST_ASSERT(1 == 0);
  //  BMP_draw_string(this->ldmap,type==0? "D'":"R^2", 10, height-50, RGB_BLACK,
  // 0);
  std::string curPath = this->path + (type == 0 ? "_D.bmp" : "_R2.bmp");
  BMP_write(this->ldmap, curPath.c_str());
  BMP_delete(this->ldmap);
  std::cout << "LD map " << (type == 0 ? "(D')" : "(R2)") << "saved to "
            << curPath << "\n";
}

void LDTest::DrawLDMapDandR2() {
  this->DrawLDMap(0);
  this->DrawLDMap(1);
}

std::string LDTest::reporthtml() {
  std::string res;
  std::string filename = get_file_name_from_full_path(this->path);
  res += "\n<h2> Linkage Disequilibrium Analysis: </h2>\n";
  res += "<h3>D'</h3>\n";
  res += "<img src=\"" + filename + "_D.bmp" + "\"" + " alt=\"LD (D') map\">\n";
  res += "<h3>R<sup>2</sup></h3>\n";
  res +=
      "<img src=\"" + filename + "_R2.bmp" + "\"" + " alt=\"LD (R2) map\">\n";
  return res;
}

} /* namespace SHEsis */
