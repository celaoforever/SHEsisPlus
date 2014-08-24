/*
 * ArrayStorage.h
 *
 *  Created on: Aug 18, 2014
 *      Author: ada
 */

#ifndef ARRAYSTORAGE_H_
#define ARRAYSTORAGE_H_
#include <boost/multi_array.hpp>
#include "Cloneable.h"
#include "utility.h"
namespace SHEsis {
class ArrayStorage: public Cloneable {
public:
	//ArrayStorage(){};
	ArrayStorage(boost::shared_ptr<int[]> sizes):sizes(sizes),dimension(sizes[0]),
		size(GeneralIndexingGetSize(sizes)),
		one_dimensional_array(new int[size+1])
	{};
	ArrayStorage(boost::shared_ptr<int[]> sizes, boost::shared_ptr<int[]> one_dimensional_array):
		sizes(sizes),dimension(sizes[0]),
			size(GeneralIndexingGetSize(sizes)),
			one_dimensional_array(one_dimensional_array)
	{};
	void set(boost::shared_ptr<int[]> indices, int value);
	void set(int index, int value);
	int get(boost::shared_ptr<int[]> indices);
	boost::shared_ptr<int[]> getArray();
	boost::shared_ptr<int[]> getSizes();
	int getDimension();
	virtual ~ArrayStorage(){};
	CLONE(ArrayStorage);
	boost::shared_ptr<int[]> sizes;
private:
	int dimension;
    int size;
	boost::shared_ptr<int[]> one_dimensional_array;
};

} /* namespace SHEsis */

#endif /* ARRAYSTORAGE_H_ */
