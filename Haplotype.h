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

struct HapTestResult{
	std::vector<boost::shared_ptr<short[]> > haplotypes;
	boost::multi_array<short,2> genotypes;
	boost::shared_ptr<int[]> CaseCount;
	boost::shared_ptr<int[]> ControlCount;
	double FisherP;
	double PearsonP;
	double ChiSquare;
	HapTestResult(int sample,int ploidy):
		genotypes(boost::extents[sample][ploidy]),FisherP(1),PearsonP(1),ChiSquare(0)
	{};
};


class Haplotype {
public:
	Haplotype(SHEsisData& data):data(data),VarNum(0),ClauseNum(0),
	occurence(boost::extents[data.getSampleNum()][data.getSnpNum()]),
	missing(boost::extents[data.getSampleNum()][data.getSnpNum()]),
	Results(data.getSampleNum(),data.getNumOfChrSet())
	//variables()
	{
		res.str("");
		sat="";
		this->statOccurence();
		for(int i=0;i<this->data.getSnpNum();i++){
				this->SnpIdx.push_back(i);
		}
	};

	Haplotype(SHEsisData& data, int Snp, std::vector<short> mask):data(data),VarNum(0),ClauseNum(0),
			occurence(boost::extents[data.getSampleNum()][Snp]),
			missing(boost::extents[data.getSampleNum()][Snp]),
			Results(data.getSampleNum(),data.getNumOfChrSet())
//			SnpIdx(Snp,0)
	{
		res.str("");
		sat="";
		this->mask=mask;
		this->statOccurenceMask();

	};

	virtual ~Haplotype(){
		for(int i=0;i<occurence.shape()[0];i++){
			for(int j=0;j<occurence.shape()[1];j++){
				this->occurence[i][j].clear();
			}
		}
		this->mask.clear();
		this->SnpIdx.clear();
		this->Results.haplotypes.clear();
	};
	SHEsisData& data;
	void statOccurence();
	void statOccurenceMask();
	void getBiallelicCoding(int ploidy,int which_genotype, int which_index, int which_allele, IndexingVariables& variables);
	void getGeneralCoding(int ploidy, int which_genotype, int which_index,IndexingVariables& variables);
	void getGeneralCodingMissing(int ploidy, int which_genotype, int which_index, IndexingVariables& variables);
	void getGeneralCodingTotalyMissing(int ploidy, int which_genotype, int which_index, IndexingVariables& variables);
	void BuildModel(IndexingVariables& variables, int number_of_explaining_haplotypes);
	void createVariables(int number_of_explaining_haplotypes,IndexingVariables& variables);
	int solve();
	void parseSolution(IndexingVariables& variables,int assumed_haplotypes);
	void associationTest();


private:
	std::stringstream res;
	std::string sat;
	int VarNum;
	int ClauseNum;
	std::vector<short> mask;
	std::vector<int> SnpIdx;

	boost::multi_array< std::vector<int>, 2> occurence;
	boost::multi_array<int, 2> missing;
	HapTestResult Results;
//	IndexingVariables variables;

};

} /* namespace SHEsis */

#endif /* HAPLOTYPE_H_ */
