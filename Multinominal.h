/*
 * Multinominal.h
 *
 *  Created on: Aug 10, 2014
 *      Author: ionadmin
 */

#ifndef MULTINOMINAL_H_
#define MULTINOMINAL_H_

#include <vector>
#include <numeric>
#include <ostream>

namespace SHEsis {
typedef std::vector<size_t> SVI;
typedef SVI::iterator SVII;
typedef SVI::const_iterator SVICI;

void view(std::ostream& ost, SVI const& v) {
  ost << "(";
  if (v.size() > 0) ost << v.at(0);

  for (size_t i = 1; i < v.size(); ++i) {
    size_t const vai(v.at(i));
    if (vai == 0) break;
    ost << ", " << vai;
  }
  ost << ")";
}

// class index calculates subscripts for class combo

class index {
  static SVI pair;
  static SVI sole;
  static size_t pair_ind(size_t const, size_t const);
  static void layer(size_t const);
  static bool compare(size_t const, size_t const);

 public:
  static size_t get(size_t const, size_t const);
  static size_t get(size_t const);
  static size_t get(SVI const&);
};

SVI index::pair(1, 1);
SVI index::sole(1, 0);

size_t index::pair_ind(size_t const rem, size_t const top) {
  return ((rem * (rem + 1)) / 2 + ((top > rem) ? rem : top));
}

void index::layer(size_t const rem) {
  if (pair.size() < pair_ind(rem, 0)) layer(rem - 1);
  pair.push_back(0);
  pair.push_back(1);
  for (size_t top = 2; top <= rem; ++top) {
    pair.push_back(pair.at(pair_ind(rem - top, top)));
    if (rem >= top) pair.back() += pair.at(pair_ind(rem, top - 1));
  }
  size_t const y(sole.back());
  sole.push_back(pair.at(pair_ind(rem - 1, rem - 1)));
  sole.back() += y;
}

bool index::compare(size_t const a, size_t const b) { return a > b; }

size_t index::get(size_t const rem, size_t const top) {
  size_t const ind(pair_ind(rem, top));
  if (ind < pair.size()) return pair.at(ind);
  layer(rem);
  return pair.at(ind);
}

size_t index::get(size_t const rem) {
  if (rem < sole.size()) return sole.at(rem);
  layer(rem);
  return sole.at(rem);
}

size_t index::get(SVI const& part) {
  SVI temp(part);
  std::sort(temp.begin(), temp.end(), &compare);
  size_t minuend(std::accumulate(temp.begin(), temp.end(), 0));
  size_t ans(get(minuend));
  for (SVICI i = temp.begin(); i != temp.end(); ++i) {
    if (*i == 0) break;
    ans += get(minuend, *i - 1);
    minuend -= *i;
  }
  return ans;
}

// the multinomial values are stored in class combo

template <typename result_type>
class combo {
  static std::vector<result_type> guts;
  static size_t tier;
  static void layer(SVI const&);
  static void layer(size_t const, size_t const, size_t const, SVI&);
  static void layer();

 public:
  static result_type get(SVI const&);
};

template <typename result_type>
std::vector<result_type> combo<result_type>::guts(1, 1);

template <typename result_type>
size_t combo<result_type>::tier(0);

template <typename result_type>
void combo<result_type>::layer(SVI const& part) {
  result_type total(0);
  SVI temp(part);
  for (SVII n = temp.begin(); n != temp.end(); ++n) {
    size_t& tan(*n);
    if (tan == 0) break;
    tan--;
    total += get(temp);
    tan++;
  }
  guts.push_back(total);
}

template <typename result_type>
void combo<result_type>::layer(size_t const rem, size_t const top,
                               size_t const pos, SVI& part) {
  if (rem == 0) {
    layer(part);
    return;
  }

  for (size_t n = 1;; ++n) {
    if (n > rem) break;
    if (n > top) break;
    part.at(pos) = n;
    if (top > n)
      layer(rem - n, n, pos + 1, part);
    else
      layer(rem - n, top, pos + 1, part);
    part.at(pos) = 0;
  }
}

template <typename result_type>
void combo<result_type>::layer() {
  ++tier;
  SVI part(tier);
  layer(tier, tier, 0, part);
}

template <typename result_type>
result_type combo<result_type>::get(SVI const& part) {
  size_t const ind(index::get(part));
  while (guts.size() <= ind) layer();
  return guts.at(ind);
}

// the next two functions are probably all that you really need

template <typename result_type>
result_type multi(SVI const& part) {
  return combo<result_type>::get(part);
}

size_t parti(size_t const rem, size_t const top) {
  return index::get(rem, (top > rem) ? rem : top);
}

template size_t multi<size_t>(SVI const& part);
}

#endif /* MULTINOMINAL_H_ */
