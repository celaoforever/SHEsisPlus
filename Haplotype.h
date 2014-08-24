/*
 * Haplotype.h
 *
 *  Created on: Aug 19, 2014
 *      Author: ada
 */

#ifndef HAPLOTYPE_H_
#define HAPLOTYPE_H_
#include "SHEsisData.h"
#include "utility.h"
#include "ArrayStorage.h"
#include "IndexingVariables.h"
namespace SHEsis {
class Haplotype {
public:
	Haplotype(SHEsisData& data):data(data),M(-1),VarNum(0),ClauseNum(0),occurence(boost::extents[data.getSampleNum()][data.getSnpNum()]){
		res.str("");
	};
	virtual ~Haplotype(){
		for(int i=0;i<data.getSampleNum();i++){
			for(int j=0;j<data.getSnpNum();j++){
				//this->occurence[i][j].clear();
			}
		}
	};
	SHEsisData& data;
	void statOccurence();
	std::string getBiallelicCoding(int ploidy,int which_genotype, int which_index, int which_allele, IndexingVariables& variables);
	std::string getGeneralCoding(int ploidy, int which_genotype, int which_index,IndexingVariables& variables);
	void BuildModel();
	IndexingVariables BuildModel(/*IndexingVariables variables_old,*/ int number_of_explaining_haplotypes);
	IndexingVariables createVariables(int);
	int statMissing(int iSample,int iSnp);
	std::vector<int> createAntiHaplotypes(int number_of_explaining_haplotypes,IndexingVariables variables,
			std::vector<int> anti_haplotypes, int number_of_known_haplotypes);
	boost::multi_array< std::vector<int>, 2> occurence;
private:
	int M;
	std::stringstream res;
	int VarNum;
	int ClauseNum;

};

} /* namespace SHEsis */

#endif /* HAPLOTYPE_H_ */
