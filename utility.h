/*
 * utility.h
 *
 *  Created on: Aug 8, 2014
 *      Author: ionadmin
 */

#ifndef UTILITY_H_
#define UTILITY_H_

#include <string>
#include <boost/assert.hpp>
#include <boost/multi_array.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#define SHEsisABS(x) x < 0 ? (-1) * x : x

typedef enum HaploMethod { EM, SAT } HapMethod;

static std::string HtmlHeader =
    "<!DOCTYPE html> \n\
	<html> \n\
	<head> \n\
	<style> \n\
	table { \n\
	    width:auto; \n\
	} \n\
	table, th, td { \n\
	    border: 1px solid black; \n\
	    border-collapse: collapse; \n\
	} \n\
	th, td { \n\
	    padding: 5px; \n\
	    text-align: left; \n\
	} \n\
	table tr:nth-child(even) { \n\
	    background-color: #eee;\n\
	}\n\
	table tr:nth-child(odd) {\n\
	   background-color:#fff;\n\
	}\n\
	table th    {\n\
	    background-color: #AAA;\n\
	    color: white;\n\
	}\n\
	</style>\n\
	</head>\n\
\n\
	<body>\n";

static std::string HtmlHeaderServer =
    "<!DOCTYPE html>\n\
		<html>\n\
		    <head>\n\
		        <link rel=\"stylesheet\" href=\"../stylesheets/metro-bootstrap.css\">\n\
			<link rel=\"stylesheet\" href=\"../stylesheets/iconFont.min.css\">\n\
		        <script src=\"../javascripts/jquery.min.js\"></script>\n\
		        <script src=\"../javascripts/jquery-ui.js\"></script>\n\
		        <script src=\"../javascripts/metro.min.js\"></script>\n\
		    </head>\n\
		    <style>\n\
			pre{\n\
			height:auto;\n\
			max-height:200px;\n\
			overflow:auto;\n\
			word-break:normal !important;\n\
			word-wrap:normal !important;\n\
			white-space:pre !important;\n\
			font-family: Arial;\n\
			padding:0;\n\
			margin:0;\n\
			-moz-tab-size:24;\n\
			-o-tab-size:24;\n\
			-webkit-tab-size:24;\n\
			-ms-tab-size:24;\n\
			tab-size:24;\n\
			}\n\
			td{\n\
			vertical-align:middle !important;\n\
			}\n\
			</style>\n\
		<body class=\"metro\">\n\
		<div class=\"navigation-bar\">\n\
		  <div class=\"navigation-bar-content container\">\n\
		    <a class=\"element\" href=\"http://analysis.bio-x.cn/SHEsisMain.htm\" title=\"Back to SHEsis main page\">\n\
		    <span class=\"icon-home\"></span>\n\
			SHEsis\n\
			<sup>2.0</sup>\n\
		     </a>\n\
		    <span class=\"element-divider\"></span>\n\
		    <a class=\"element\" href=\"../Help.html\">\n\
		      Help\n\
		    </a>\n\
		    <span class=\"element-divider\"></span>\n\
		    <a class=\"element\" href=\"../CiteUs.html\">\n\
		      Cite Us \n\
		    </a>   \n\
		    <span class=\"element-divider\"></span>\n\
		    <a class=\"element\" href=\"../ContactUs.html\">\n\
		      Contact Us \n\
		    </a> \n\
		    <span class=\"element-divider\"></span>\n\
		    <a class=\"element place-right\" href=\"https://github.com/celaoforever/SHEsis\" title=\"Follow us on github\">\n\
		      Follow us on github\n\
		      <span class=\"icon-github-2\"> \n\
		    </a>\n\
		 </div>\n\
		</div>\n\
		<div class=\"container\"  style=\"padding:20px 50px\">\n\
		<h1>\n\
		 <a href=\"/\">\n\
		  <i class=\"icon-arrow-left-3 fg-darker smaller\"></i>\n\
		 </a>\n\
		  Results\n\
		 <small class=\"on-right\"></small>\n\
		</h1>\n\
";

template <class BidIt>
inline bool next_combination(BidIt n_begin, BidIt n_end, BidIt r_begin,
                             BidIt r_end) {

  bool boolmarked = false;
  BidIt r_marked;

  BidIt n_it1 = n_end;
  --n_it1;

  BidIt tmp_r_end = r_end;
  --tmp_r_end;

  for (BidIt r_it1 = tmp_r_end; r_it1 != r_begin || r_it1 == r_begin;
       --r_it1, --n_it1) {
    if (*r_it1 == *n_it1) {
      if (r_it1 != r_begin)  // to ensure not at the start of r sequence
      {
        boolmarked = true;
        r_marked = (--r_it1);
        ++r_it1;  // add it back again
        continue;
      } else  // it means it is at the start the sequence, so return false
        return false;
    } else  // if(*r_it1!=*n_it1 )
    {
      // marked code
      if (boolmarked == true) {
        // for loop to find which marked is in the first sequence
        BidIt n_marked;  // mark in first sequence
        for (BidIt n_it2 = n_begin; n_it2 != n_end; ++n_it2)
          if (*r_marked == *n_it2) {
            n_marked = n_it2;
            break;
          }

        BidIt n_it3 = ++n_marked;
        for (BidIt r_it2 = r_marked; r_it2 != r_end; ++r_it2, ++n_it3) {
          *r_it2 = *n_it3;
        }
        return true;
      }
      for (BidIt n_it4 = n_begin; n_it4 != n_end; ++n_it4)
        if (*r_it1 == *n_it4) {
          *r_it1 = *(++n_it4);
          return true;
        }
    }
  }

  return true;  // will never reach here
}

template <class T>
class array1D {
 public:
  array1D() {};
  array1D(int msize) : isStatic(false), msize(msize) {
    this->a = new T[msize];
    std::memset(a, 0, msize * sizeof(T));
  }
  array1D(T* a, int msize) : a(a), msize(msize), isStatic(true) {
    this->a = new T[msize];
    // for(int )
  };
  ~array1D() {
    delete[] a;
    a = 0;
  }
  int size() {
    return msize;
  };
  T& operator[](const int idx) {
    BOOST_ASSERT(idx < msize);
    return a[idx];
  }

 private:
  T* a;
  int msize;
  bool isStatic;
};
template <typename T>
std::string convert2string(T v) {
  std::stringstream ss;
  if (v == -999) return "NA";
  if (v == 0) return "0";
  if (v > 0.001 || v < -0.001) {
    int tmp = v * 1000;
    v = (T)tmp / 1000.0;
    ss << v;
  } else {
    ss << std::setprecision(2) << std::scientific << v;
  }

  return ss.str();
}

struct MultiComp {
  double p;
  int idx;
  bool operator<(const MultiComp& p2) const { return (p < p2.p); }
};

int getRank(double p, std::vector<double> v);
std::string get_file_name_from_full_path(const std::string& file_path);
std::string ToBinaryString(int i);
boost::shared_ptr<int[]> SetSharedPtr(int Num, ...);
std::string int2str(int n);
void error(std::string msg);
int GeneralIndexingGetIndex(boost::shared_ptr<int[]> sizes,
                            boost::shared_ptr<int[]> indices);
int GeneralIndexingGetSize(boost::shared_ptr<int[]> sizes);
boost::shared_ptr<int[]> toBooleanInt(int n, int i);
void HolmCorrection(std::vector<MultiComp>& p, std::vector<double>& adjusted);
void SidakSSCorrection(std::vector<MultiComp>& p,
                       std::vector<double>& adjusted);
void SidakSDCorrection(std::vector<MultiComp>& p,
                       std::vector<double>& adjusted);
void BHCorrection(std::vector<MultiComp>& p, std::vector<double>& adjusted);
void BYCorrection(std::vector<MultiComp>& p, std::vector<double>& adjusted);

template <typename T>
// tab is stored in column-major format
void PearsonChiSquareTest(T* tab, int row, int col, double& chi, double& p) {
  boost::multi_array<T, 1> RowSum(boost::extents[row]);
  std::fill(RowSum.begin(), RowSum.end(), 0);

  boost::multi_array<T, 1> ColSum(boost::extents[col]);
  std::fill(ColSum.begin(), ColSum.end(), 0);

  int degreeOfFreedom = (row - 1) * (col - 1);
  if (degreeOfFreedom <= 0) {
    chi = -999;
    p = -999;
    return;
  }
  T total = 0;
  chi = 0;
  for (int i = 0; i < RowSum.shape()[0]; i++)
    for (int j = 0; j < ColSum.shape()[0]; j++) {
      T observe = tab[j * row + i];
      RowSum[i] += observe;
      ColSum[j] += observe;
      total += observe;
    }
  if (0 == total) {
    chi = 0;
    p = 1;
    return;
  }

  for (int i = 0; i < RowSum.shape()[0]; i++)
    for (int j = 0; j < ColSum.shape()[0]; j++) {
      T observe = tab[j * row + i];
      T expected = RowSum[i] * ColSum[j] / total;
      if (0 == expected) continue;
      chi += (observe - expected) * (observe - expected) / expected;
    }
  try {
    boost::math::chi_squared dist(degreeOfFreedom);
    p = boost::math::cdf(boost::math::complement(dist, chi));
  }
  catch (...) {
    p = -999;
  }
}

template void PearsonChiSquareTest<double>(double* tab, int row, int col,
                                           double& chi, double& p);
#endif /* UTILITY_H_ */
