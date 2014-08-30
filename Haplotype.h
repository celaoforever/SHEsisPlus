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
	Haplotype(SHEsisData& data):data(data),VarNum(0),ClauseNum(0),
	occurence(boost::extents[data.getSampleNum()][data.getSnpNum()]),
	missing(boost::extents[data.getSampleNum()][data.getSnpNum()])
	{
		res.str("");
		this->statOccurence();
	};

	Haplotype(SHEsisData& data, int Snp, std::vector<short> mask):data(data),VarNum(0),ClauseNum(0),
			occurence(boost::extents[data.getSampleNum()][Snp]),
			missing(boost::extents[data.getSampleNum()][Snp])
	{
		res.str("");
		this->mask=mask;
		this->statOccurenceMask();
	};

	virtual ~Haplotype(){
		for(int i=0;i<data.getSampleNum();i++){
			for(int j=0;j<data.getSnpNum();j++){
				this->occurence[i][j].clear();
			}
		}
		this->mask.clear();
	};
	SHEsisData& data;
	void statOccurence();
	void statOccurenceMask();
	void getBiallelicCoding(int ploidy,int which_genotype, int which_index, int which_allele, IndexingVariables& variables);
	void getGeneralCoding(int ploidy, int which_genotype, int which_index,IndexingVariables& variables);
	void getGeneralCodingMissing(int ploidy, int which_genotype, int which_index, IndexingVariables& variables);
	void getGeneralCodingTotalyMissing(int ploidy, int which_genotype, int which_index, IndexingVariables& variables);
	void BuildModel();
	void BuildModel(/*IndexingVariables variables_old,*/ int number_of_explaining_haplotypes);
	void createVariables(int number_of_explaining_haplotypes,IndexingVariables& variables);
	int statMissing(int iSample,int iSnp);
	boost::multi_array< std::vector<int>, 2> occurence;
	boost::multi_array<int, 2> missing;
private:
	std::stringstream res;
	int VarNum;
	int ClauseNum;
	std::vector<short> mask;

};

} /* namespace SHEsis */

#endif /* HAPLOTYPE_H_ */
