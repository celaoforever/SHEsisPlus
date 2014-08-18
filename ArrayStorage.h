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
template <class T>
class ArrayStorage: public Cloneable {
public:
	ArrayStorage(boost::multi_array<int,1>& sizes):sizes(sizes),dimension(sizes.shape()[0]),
		size(GeneralIndexingGetSize(sizes)),one_dimensional_array(boost::extents[size])
	{};
	ArrayStorage(boost::multi_array<int,1> sizes, boost::multi_array<T,1> one_dimensional_array):
		sizes(sizes),dimension(sizes.shape()[0]),
			size(GeneralIndexingGetSize(sizes)),one_dimensional_array(one_dimensional_array)
	{};

	void set(boost::multi_array<int,1>& indices, int value){
        int index = GeneralIndexingGetIndex(sizes, indices);
        BOOST_ASSERT(index < size);
        one_dimensional_array[index] = value;
	};

	void set(int index, int value){
		one_dimensional_array[index] = value;
	};

	T get(boost::multi_array<int,1>& indices){
        int index = GeneralIndexingGetIndex(sizes, indices);
        BOOST_ASSERT(index < size);
        return one_dimensional_array[index];
	};
	boost::multi_array<T,1>& getArray(){
		return this->one_dimensional_array;
	}
	boost::multi_array<int,1>& getSizes(){
		return this->sizes;
	}
	int getDimension(){
		return this->dimension;
	};
	virtual ~ArrayStorage();
	CLONE(ArrayStorage);
private:
	int dimension;
	boost::multi_array<T,1> one_dimensional_array;
    boost::multi_array<int,1> sizes;
    int size;
};

} /* namespace SHEsis */

#endif /* ARRAYSTORAGE_H_ */
