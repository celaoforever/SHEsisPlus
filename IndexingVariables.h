/*
 * IndexingVariables.h
 *
 *  Created on: Aug 18, 2014
 *      Author: ada
 */

#ifndef INDEXINGVARIABLES_H_
#define INDEXINGVARIABLES_H_
#include "utility.h"
#include "ArrayStorage.h"
#include <boost/unordered_map.hpp>
namespace SHEsis{
template <class T, class S, class R>
class IndexingVariables : public Cloneable{
public:

	CLONE(IndexingVariables);
	IndexingVariables(){
		this->reset();
		this->offset=1;
	}
	~IndexingVariables(){
		this->hm.clear();
		this->ILPsolution.clear();
		this->id.clear();
		this->enumeration.clear();
	}
	boost::unordered_map<int,int> id;
	boost::unordered_map<int,int> enumeration;
	void reset(){
	   this->id.clear();
	   this->enumeration.clear();
	}
	void add(T key, boost::multi_array<int,1> sizes)
	{
	    int size = GeneralIndexingGetSize(sizes);
	    boost::multi_array<int,1> variables(boost::extents[size]);
	    for(int i = 0; i < size; i++)
	        variables[i] = offset + i;
	    offset += size;
	    ArrayStorage as(sizes, variables);
	    hm.put(key, as);
	};

	int getEnumeration(int ident){
		BOOST_ASSERT(ident>=1);
		return this->id[ident];
	};

	int getEnumeration(T key, boost::multi_array<int ,1> indices){
		int ident=SHEsisABS((this->hm[key]).get(indices));
		int size=0;
		if(id.end() == id.find(ident)){
			size = id.size()+1;
			id[ident]=size;
			enumeration[size]=ident;
			return size;
		}else{
			return id[ident];
		}
	};

	std::string getEnumerationILP( T key, boost::multi_array<int ,1> indices){
		int ident=SHEsisABS((this->hm[key]).get(indices));
		int size=0;
		std::stringstream res;
		if(id.end() == id.find(ident)){
			size=id.size()+1;
			id[ident]=size;
			enumeration[size]=ident;
			res<<"ilp."<<size;
			return res.str();
		}else
		{
			res<<"ilp."<<id[ident];
			return res.str();
		}
	};

	int getId(int enumerat){
		BOOST_ASSERT(enumerat>=1);
		return enumeration[enumerat];
	}

	int getEvalutatedId(T key, boost::multi_array<int, 1> indices){
		return this->hm[key].get(indices);
	};

	double getEvalutatedIdILP(T key, boost::multi_array<int,1> indices){
		int ident=this->hm[key].get(indices);
		if(id.end() != id.find(ident)){
			std::stringstream k;
			k<<"ilp."<<id[ident];
			return ILPsolution[k.str()];
		}else
		{
			return NULL;
		}
	};

	void set( T key, boost::multi_array<S,1> variables, boost::multi_array<int,1> sizes ){
		ArrayStorage<S> as(sizes, variables);
		hm[key]=as;
	};

	void set(boost::unordered_map<T,ArrayStorage>& as, int offset){
		hm.clear();
		hm=as;
		this->offset=offset;
	}

	boost::multi_array<int, 1>& getObjectSizes(T key){
		if(hm.end() != hm.find(key))
			return hm[key].getSizes();
		return NULL;
	};

	int size()
	{
		return hm.size();
	};

	int numberOfVariables()
	{
		return this->offset-1;
	}

	void setParities(boost::multi_array<int,1> parities, boost::multi_array<R, 1> outsider){
		//TO-DO
	}

	void setValueILP(std::string key, double value){
		ILPsolution[key]=value;
	};


	void show(){};
	int test(){return id.size();};

	int getEnumerationSize(){
		return enumeration.size();
	};

	bool containsKey(T key){
		return (hm.end() != hm.find(key));
	};

	int getIdSize(){
		return id.size();
	};

private:
	boost::unordered_map<T,ArrayStorage> hm;
	boost::unordered_map<std::string,double> ILPsolution;
	int offset;

};

}
#endif /* INDEXINGVARIABLES_H_ */
