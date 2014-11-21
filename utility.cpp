/*
 * utility.cpp
 *
 *  Created on: Aug 8, 2014
 *      Author: ionadmin
 */
#include "utility.h"
#include <boost/array.hpp>
#include <stdarg.h>
#include <string>
#include <algorithm>

void HolmCorrection(std::vector<MultiComp>& p, std::vector<double>& adjusted) {
  // p is sorted
  int count = p.size();
  adjusted.resize(count, 0);
  if (p[0].p < 0)
    adjusted[p[0].idx] = -999;
  else
    adjusted[p[0].idx] = (p[0].p * count > 1 ? 1 : p[0].p * count);
  for (int i = 1; i < count; i++) {
    if (p[i].p < 0) {
      adjusted[p[i].idx] = -999;
      continue;
    };
    double x = (count - i) * p[i].p < 1 ? (count - i) * p[i].p : 1;
    adjusted[p[i].idx] =
        adjusted[p[i - 1].idx] > x ? adjusted[p[i - 1].idx] : x;
  };
}

void SidakSSCorrection(std::vector<MultiComp>& p,
                       std::vector<double>& adjusted) {
  int count = p.size();
  adjusted.resize(count, 0);
  for (int i = 0; i < count; i++) {
    if (p[i].p < 0) {
      adjusted[p[i].idx] = -999;
      continue;
    };
    adjusted[p[i].idx] = 1 - pow(1 - p[i].p, count);
  }
}

void SidakSDCorrection(std::vector<MultiComp>& p,
                       std::vector<double>& adjusted) {
  int count = p.size();
  adjusted.resize(count, 0);
  if (p[0].p < 0)
    adjusted[p[0].idx] = -999;
  else
    adjusted[p[0].idx] = 1 - pow(1 - p[0].p, count);
  for (int i = 1; i < count; i++) {
    if (p[i].p < 0) {
      adjusted[p[i].idx] = -999;
      continue;
    };
    double x = 1 - pow(1 - p[i].p, count - i);
    adjusted[p[i].idx] =
        adjusted[p[i - 1].idx] > x ? adjusted[p[i - 1].idx] : x;
  }
}

void BHCorrection(std::vector<MultiComp>& p, std::vector<double>& adjusted) {
  int count = p.size();
  adjusted.resize(count, 0);
  if (p[count - 1].p < 0)
    adjusted[p[count - 1].idx] = -999;
  else
    adjusted[p[count - 1].idx] = p[count - 1].p;
  for (int i = count - 2; i >= 0; i--) {
    if (p[i].p < 0) {
      adjusted[p[i].idx] = -999;
      continue;
    };
    double x = ((double)count) / ((double)(i + 1)) * p[i].p < 1
                   ? ((double)count / ((double)(i + 1))) * p[i].p
                   : 1;
    adjusted[p[i].idx] =
        adjusted[p[i + 1].idx] < x ? adjusted[p[i + 1].idx] : x;
  }
}

void BYCorrection(std::vector<MultiComp>& p, std::vector<double>& adjusted) {
  int count = p.size();
  adjusted.resize(count, 0);
  double a = 0;
  for (double i = 1; i <= count; i++) a += 1 / i;
  if (p[count - 1].p < 0)
    adjusted[p[count - 1].idx] = -999;
  else
    adjusted[p[count - 1].idx] =
        a * p[count - 1].p < 1 ? a * p[count - 1].p : 1;

  for (int i = count - 2; i >= 0; i--) {
    if (p[i].p < 0) {
      adjusted[p[i].idx] = -999;
      continue;
    };
    double x = (((double)count * a) / (double(i + 1))) * p[i].p < 1
                   ? (((double)count * a) / ((double)(i + 1))) * p[i].p
                   : 1;
    adjusted[p[i].idx] =
        adjusted[p[i + 1].idx] < x ? adjusted[p[i + 1].idx] : x;
  }
};

int getRank(double p, std::vector<double> v) {
  if (p < v[0]) return 0;
  for (int i = 0; i < v.size() - 1; i++) {
    if (p >= v[i] && p <= v[i + 1]) return i;
  }
  return v.size();
}

boost::shared_ptr<int[]> SetSharedPtr(int Num, ...) {
  boost::shared_ptr<int[]> sp(new int[Num + 1]);
  sp[0] = Num;
  va_list args;
  va_start(args, Num);
  int idx = 1;
  while (Num > 0) {
    sp[idx++] = va_arg(args, int);
    Num--;
  }
  return sp;
}

std::string get_file_name_from_full_path(const std::string& file_path) {
  std::string file_name;

  std::string::const_reverse_iterator it =
      std::find(file_path.rbegin(), file_path.rend(), '/');
  if (it != file_path.rend()) {
    file_name.assign(file_path.rbegin(), it);
    std::reverse(file_name.begin(), file_name.end());
    return file_name;
  } else
    return file_path;
}

char digits[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b',
                 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
                 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'};

std::string ToUnsignedString(int i, int shift) {
  char buf[32];
  char* pBuf = buf;
  int charPos = 32;
  int radix = 1 << shift;
  int mask = radix - 1;
  do {
    pBuf[--charPos] = digits[i & mask];
    i = i >> shift;
  } while (i != 0);

  std::string str;
  int strLen = 32 - charPos;
  pBuf = pBuf + charPos;
  while (strLen) {
    str.push_back(*pBuf);
    pBuf++;
    strLen--;
  }

  return str;
}

std::string ToBinaryString(int i) { return ToUnsignedString(i, 1); }
std::string int2str(int n) {
  std::ostringstream s2(std::stringstream::out);
  s2 << n;
  return s2.str();
}

boost::shared_ptr<int[]> toBooleanInt(int n, int i) {

  // boost::shared_ptr<int[]> spb(new int[n+1]);
  boost::shared_ptr<int[]> spb(new int[n]);
  for (int i = 0; i < n; i++) {
    spb[i] = 0;
  }
  boost::shared_ptr<int[]> spb_(new int[n]);
  for (int i = 0; i < n; i++) {
    spb_[i] = 0;
  }
  std::string s = ToBinaryString(i);
  for (int j = n - s.length(); j < n; j++) {
    char c = s[j - n + s.length()];
    if (c == '1') spb_[j] = 1;
  }
  int counter = n - 1;
  for (int j = 0; j < n; j++) {
    spb[j] = spb_[counter];
    counter--;
  }
  return spb;
}

void error(std::string msg) {
  //	std::cout<<msg<<"\n";
  throw std::runtime_error("error in calculating fisher's exact test");
}
/* translated from clase GenralIndexing*/
int GeneralIndexingGetIndex(boost::shared_ptr<int[]> sizes,
                            boost::shared_ptr<int[]> indices) {
  int index = 0;
  for (int i = 1; i <= indices[0]; i++) {
    int p = 1;
    for (int j = i + 1; j <= sizes[0]; j++) p *= sizes[j];

    p *= indices[i];
    index += p;
  }
  return ++index;
};

int GeneralIndexingGetSize(boost::shared_ptr<int[]> sizes) {
  int p = 1;
  for (int i = 1; i <= sizes[0]; i++) p *= sizes[i];
  return p;
};

// template<typename T>
//
// double CalChiSquareFromContigencyTable(boost::multi_array<T, 2> &array){
//	boost::multi_array<T,1> RowSum;
//	RowSum(boost::extents[array.shape()[0]]);
//	std::fill(RowSum.begin(),RowSum.end(),0);
//
//	boost::multi_array<T,1> ColSum;
//	ColSum(boost::extents[array.shape()[1]]);
//	std::fill(ColSum.begin(),ColSum.end(),0);
//
//	T total=0;
//	double chi=0;
//	for(int i=0;i<RowSum.shape()[0];i++)
//		for(int j=0; j<ColSum.shape()[0];j++)
//		{
//			T observe=array[i][j];
//			RowSum[i]+=observe;
//			ColSum[j]+=observe;
//			total+=observe;
//		}
//	if(0 == total)
//		return 0;
//
//	for(int i=0;i<RowSum.shape()[0];i++)
//		for(int j=0; j<ColSum.shape()[0];j++)
//		{
//			T observe=array[i][j];
//			T expected=RowSum[i]*ColSum[j]/total;
//			if(0 == expected )
//				continue;
//			chi+=(observe-expected)*(observe-expected)/expected;
//		}
//
//	return chi;
//}
//
////template double
/// CalChiSquareFromContigencyTable<double>(boost::multi_array<double, 2>
/// array);
//
//
////Copied from fisher.c - Coded by Matthew Belmonte
// double fisher(int m, int n, double x)
//{
//	int a, b, i, j;
//	double w, y, z, zk, d, p;
//	a = 2*(m/2)-m+2;
//	b = 2*(n/2)-n+2;
//	w = (x*m)/n;
//	z = 1.0/(1.0+w);
//	if(a == 1)
//	{
//		if(b == 1)
//		{
//			p = sqrt(w);
//			y = 0.3183098862;
//			d = y*z/p;
//			p = 2.0*y*atan(p);
//		}
//		else
//		{
//			p = sqrt(w*z);
//			d = 0.5*p*z/w;
//		}
//	}
//	else if(b == 1)
//	{
//		p = sqrt(z);
//		d = 0.5*z*p;
//		p = 1.0-p;
//	}
//	else
//	{
//		d = z*z;
//		p = w*z;
//	}
//	y = 2.0*w/z;
//	if(a == 1)
//		for(j = b+2; j <= n; j += 2)
//		{
//			d *= (1.0+1.0/(j-2))*z;
//			p += d*y/(j-1);
//		}
//	else
//	{
//		zk = pow(z, (double)((n-1)/2));
//		d *= (zk*n)/b;
//		p = p*zk+w*z*(zk-1.0)/(z-1.0);
//	}
//	y = w*z;
//	z = 2.0/z;
//	b = n-2;
//	for(i = a+2; i <= m; i += 2)
//	{
//		j = i+b;
//		d *= (y*j)/(i-2);
//		p -= z*d/j;
//	}
//	return(p<0.0? 0.0: p>1.0? 1.0: p);
//}
//
// double getFisherP(double chi2, int df)
//{
//	if(fabs(chi2)<1e-10) return 1.0;
//	else return 1.0-fisher(df, 5000, chi2/df);
//}
//
//// upper one sided tail probability of the normal distribution for a given
/// normal deviate
// double probnorm(double x)
//{
//	double z,t,p,xa;
//	xa=fabs(x);
//	if(xa>12.0) p=0.0;
//	else
//	{
//		z=0.39894228*exp(-0.5*xa*xa);
//		t=1.0/(1.0+0.33267*xa);
//		p=z*t*(0.4361836+t*(0.937298*t-0.1201676));
//	}
//	if(x>=0.0) return p;
//	else return(1.0-p);
//}
//
//
//// Pearson's chi2prob
// double getPearsonP(double x2,int ndf)
//{
//	double z,p,x,sum,re,ch,chp;
//	short i,n1,n2;
//
//	if(x2<=0.0&&ndf<=0)
//		return 1.0;
//	/*ndf=1*/
//	if(ndf==1)
//	{
//		x=sqrt(x2);return(2.0*probnorm(x));
//	}
//	/*ndf>1*/
//	if(x2>169.0)
//	{                   /*formula 26.4.14, p.941*/
//		ch=2.0/(9.0*ndf);
//		x=(exp(log(x2/ndf)/3.0)-1.0+ch)/sqrt(ch);
//		return(probnorm(x));
//	}
//	/*ndf=2*/
//	if(ndf==2)
//		return(exp(-0.5*x2));
//	/*ndf>2*/
//	n1=(ndf-1)/2;
//	n2=(ndf-2)/2;
//	if(n1==n2) /*ndf is even and >2*/
//	{
//		sum=1.0;
//		re=0.5*x2;
//		ch=1.0;
//		for(i=1; i<=n2; i++)
//		{
//			ch=ch*re/i;
//			sum+=ch;
//		}
//		return(exp(-re)*sum);
//	}
//	else
//		/*ndf is odd and >1*/
//	{
//		ch=sqrt(x2);
//		z=0.39894228*exp(-0.5*x2);
//		p=probnorm(ch);
//		if(ndf==3) return(2.0*(p+z*ch)); /*ndf=3*/
//		else     /*ndf odd and >3*/
//		{
//			chp=ch;
//			re=1.0;
//			for(i=2; i<=n1; i++)
//			{
//				re+=2.0;
//				chp=chp*x2/re;
//				ch+=chp;
//			}
//			return(2.0*(p+z*ch));
//		}
//	}
//}
