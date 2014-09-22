/*
 * ArrayStorage.cpp
 *
 *  Created on: Aug 18, 2014
 *      Author: ada
 */

#include "ArrayStorage.h"
namespace SHEsis {

void ArrayStorage::set(boost::shared_ptr<int[]> indices, int value) {
  int index = GeneralIndexingGetIndex(sizes, indices);
  BOOST_ASSERT(index <= one_dimensional_array[0]);
  one_dimensional_array[index] = value;
};

void ArrayStorage::set(int index, int value) {
  one_dimensional_array[index] = value;
};

int ArrayStorage::get(boost::shared_ptr<int[]> indices) {
  int index = GeneralIndexingGetIndex(sizes, indices);
  BOOST_ASSERT(index <= one_dimensional_array[0]);
  return one_dimensional_array[index];
};
boost::shared_ptr<int[]> ArrayStorage::getArray() {
  return this->one_dimensional_array;
}
boost::shared_ptr<int[]> ArrayStorage::getSizes() { return this->sizes; }
int ArrayStorage::getDimension() {
  return this->dimension;
};

} /* namespace SHEsis */
