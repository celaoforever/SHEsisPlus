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
	void GenerateSNPCombination(int snpnum,std::vector<std::vector<int> >& ret);
	void GenerateSNPCombination(int snpnum,std::vector<int> snpidx,std::vector<std::vector<int> >& ret);
	void GenerateGenotypeCombination(std::vector<int>& Snp,std::vector<std::vector<std::string> > & ret );
	bool genotypeEqual(std::string geno,int sample,int snp);
	double getEntropy(std::vector<int> samples, std::vector<int> snps);
	double getInformationInteraction(std::vector<int> samples,std::vector<int> snps);
};

} /* namespace SHEsis */

#endif /* GENEINTERACTION_H_ */
