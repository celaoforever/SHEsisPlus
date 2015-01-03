/*
 * GeneInteraction.h
 *
 *  Created on: Dec 30, 2014
 *      Author: ionadmin
 */

#ifndef GENEINTERACTION_H_
#define GENEINTERACTION_H_
#include "SHEsisData.h"
#include "boost/shared_ptr.hpp"
#include "utility.h"
template bool next_combination(std::vector<int>::iterator n_begin,
                               std::vector<int>::iterator n_end,
                               std::vector<int>::iterator r_begin,
                               std::vector<int>::iterator r_end);
namespace SHEsis {
typedef enum{
	EQUAL=0,
	NOT_EQUAL=1,
	MISS=2
}genotypeRetCode;


class GeneInteraction {
public:
	boost::shared_ptr<SHEsisData> data;
	GeneInteraction(boost::shared_ptr<SHEsisData> data);
	virtual ~GeneInteraction();
	void setlb(int l){this->lowbound=l;}
	void setub(int u){this->upperbound=u;}
	virtual void CalGeneInteraction(){};
protected:
	int lowbound;
	int upperbound;
	boost::multi_array<std::string,2> mGenotypeStr;
	genotypeRetCode genotypeEqual(std::string geno,int sample,int snp);
	virtual double GetSingleSNPEntroyp(int snp,int group){};
	virtual double GetMutualEntropy(std::vector<int> Snp){};
	virtual void  GenerateGenotypeCombination(std::vector<int> Snp,std::vector<std::vector<std::string> >& ret ){};
	virtual double GenotypeCombinationCount(std::vector<int> Snp,std::vector<int> Genotype,int group){};

};

} /* namespace SHEsis */

#endif /* GENEINTERACTION_H_ */
