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
namespace SHEsis {
class IndexingVariables : public Cloneable {
 public:
  CLONE(IndexingVariables);
  IndexingVariables() { this->reset(); }
  ~IndexingVariables() {
    this->hm.clear();
    this->ILPsolution.clear();
    this->id.clear();
    this->enumeration.clear();
  }
  void printVar() {
    boost::unordered_map<std::string,
                         boost::shared_ptr<ArrayStorage> >::iterator iter;
    for (iter = this->hm.begin(); iter != this->hm.end(); iter++) {
      std::cout << iter->first << "\n";
      iter->second->printSizeAndVar();
      std::cout << "\n";
    }
  }
  boost::unordered_map<int, int> id;
  boost::unordered_map<int, int> enumeration;
  void reset();
  void add(std::string key, boost::shared_ptr<int[]> sizes);
  int getEnumeration(int ident);
  int getEnumeration(std::string key, boost::shared_ptr<int[]> indices);
  std::string getEnumerationILP(std::string key,
                                boost::shared_ptr<int[]> indices);
  int getId(int enumerat);
  int getEvalutatedId(std::string key, boost::shared_ptr<int[]> indices);
  double getEvalutatedIdILP(std::string key, boost::shared_ptr<int[]> indices);
  void set(std::string key, boost::shared_ptr<int[]> variables,
           boost::shared_ptr<int[]> sizes);
  void set(
      boost::unordered_map<std::string, boost::shared_ptr<ArrayStorage> > as,
      int offset);
  boost::shared_ptr<int[]> getObjectSizes(std::string key);
  int size();
  void setParities(boost::shared_ptr<int[]> parities);
  void setValueILP(std::string key, double value);
  void show();
  int test();
  int getEnumerationSize();
  bool containsKey(std::string key);
  int getIdSize();
  void printHmKey();
  int getVarnum() {
    return this->varnum;
  };
  int numberOfVariables() {
    return offset - 1;
  };

 private:
  boost::unordered_map<std::string, boost::shared_ptr<ArrayStorage> > hm;
  boost::unordered_map<std::string, double> ILPsolution;
  int offset;
  int varnum;
};
}
#endif /* INDEXINGVARIABLES_H_ */
