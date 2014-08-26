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
	void printVar(){
		boost::unordered_map<std::string,boost::shared_ptr<ArrayStorage> > ::iterator iter;
		for(iter=this->hm.begin();iter!=this->hm.end();iter++){
			std::cout<<iter->first<<"\n";
			iter->second->printSizeAndVar();
			std::cout<<"\n";
		}
	}
	boost::unordered_map<int,int> id;
	boost::unordered_map<int,int> enumeration;
	void reset();
	void add(std::string key, boost::shared_ptr<int[]> sizes);
	int getEnumeration(int ident);
	int getEnumeration(std::string key,boost::shared_ptr<int[]> indices);
	std::string getEnumerationILP( std::string key, boost::shared_ptr<int[]> indices);
	int getId(int enumerat);
	int getEvalutatedId(std::string key, boost::shared_ptr<int[]>indices);
	double getEvalutatedIdILP(std::string key, boost::shared_ptr<int[]> indices);
	void set( std::string key, boost::shared_ptr<int[]> variables, boost::shared_ptr<int[]> sizes );
	void set(boost::unordered_map<std::string,boost::shared_ptr<ArrayStorage> > as, int offset);
	boost::shared_ptr<int[]> getObjectSizes(std::string key);
	int size();

	int numberOfVariables();
	void setParities(boost::multi_array<int,1> parities, boost::multi_array<int, 1> outsider);
	void setValueILP(std::string key, double value);
	void show();
	int test();
	int getEnumerationSize();
	bool containsKey(std::string key);
	int getIdSize();
	void printHmKey();
private:
	boost::unordered_map<std::string,boost::shared_ptr<ArrayStorage> > hm;
	boost::unordered_map<std::string,double> ILPsolution;
	int offset;

};

}
#endif /* INDEXINGVARIABLES_H_ */
