/*
 * Haplotype.h
 *
 *  Created on: Aug 19, 2014
 *      Author: ada
 */

#ifndef HAPLOTYPE_H_
#define HAPLOTYPE_H_
#include "HaplotypeBase.h"
#include "utility.h"
#include "ArrayStorage.h"
#include "IndexingVariables.h"
namespace SHEsis {

class Haplotype : public HaplotypeBase {
 public:
  Haplotype(boost::shared_ptr<SHEsisData> data);
  Haplotype(boost::shared_ptr<SHEsisData> data, int Snp,
            std::vector<short> mask);
  virtual ~Haplotype();
  virtual void startHaplotypeAnalysis();

 private:
  std::stringstream res;
  std::string sat;
  int VarNum;
  int ClauseNum;
  boost::multi_array<std::vector<int>, 2> occurence;
  boost::multi_array<int, 2> missing;
  boost::shared_ptr<IndexingVariables> variables;
  void statOccurence();
  void statOccurenceMask();
  void getBiallelicCoding(int ploidy, int which_genotype, int which_index,
                          int which_allele);
  void getGeneralCoding(int ploidy, int which_genotype, int which_index);
  void getGeneralCodingMissing(int ploidy, int which_genotype, int which_index);
  void getGeneralCodingTotalyMissing(int ploidy, int which_genotype,
                                     int which_index);
  void BuildModel(int number_of_explaining_haplotypes);
  void createVariables(int number_of_explaining_haplotypes);
  int solve();
  void parseSolution(int assumed_haplotypes);
};

} /* namespace SHEsis */

#endif /* HAPLOTYPE_H_ */
