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
  int tmp = v * 1000;
  v = (T)tmp / 1000.0;
  if (v != -1)
    ss << v;
  else
    ss << "NA";
  return ss.str();
}
std::string get_file_name_from_full_path(const std::string& file_path);
std::string ToBinaryString(int i);
boost::shared_ptr<int[]> SetSharedPtr(int Num, ...);
std::string int2str(int n);
void error(std::string msg);
int GeneralIndexingGetIndex(boost::shared_ptr<int[]> sizes,
                            boost::shared_ptr<int[]> indices);
int GeneralIndexingGetSize(boost::shared_ptr<int[]> sizes);
boost::shared_ptr<int[]> toBooleanInt(int n, int i);

template <typename T>
// tab is stored in column-major format
void PearsonChiSquareTest(T* tab, int row, int col, double& chi, double& p) {
  boost::multi_array<T, 1> RowSum(boost::extents[row]);
  std::fill(RowSum.begin(), RowSum.end(), 0);

  boost::multi_array<T, 1> ColSum(boost::extents[col]);
  std::fill(ColSum.begin(), ColSum.end(), 0);

  int degreeOfFreedom = (row - 1) * (col - 1);
  if(degreeOfFreedom<=0){
	  	 chi=-1;
	  	 p=-1;
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
  boost::math::chi_squared dist(degreeOfFreedom);
  p = boost::math::cdf(boost::math::complement(dist, chi));
}

template void PearsonChiSquareTest<double>(double* tab, int row, int col,
                                           double& chi, double& p);
#endif /* UTILITY_H_ */
