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
#include "CreatHtmlTable.h"
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
	void setlb(int l){this->lowbound=l<2?2:l;}
	void setub(int u){if(this->data && u>this->data->getSnpNum()){this->upperbound=this->data->getSnpNum();}else{this->upperbound=u;}}
	void setAdjust(bool s){this->adjust=s;};
	void setPermutation(int s){this->permutation=s;};
	virtual void setBinNum(int s){};
	virtual void setMinBin(int s){};
	virtual void setSamplePerBin(int s){};
	virtual void CalGeneInteraction(){};
	virtual std::string reporttxt(){};
	virtual std::string reporthtml(){};
protected:
	int lowbound;
	int upperbound;
	bool adjust;
	int permutation;
	boost::multi_array<std::string,2> mGenotypeStr;
	double getP(std::vector<double> permutated,double origin);
	void GenerateSNPCombination(int snpnum,std::vector<std::vector<int> >& ret);
	void GenerateSNPCombination(int snpnum,std::vector<int> snpidx,std::vector<std::vector<int> >& ret);
	void GenerateGenotypeCombination(std::vector<int>& Snp,std::vector<std::vector<std::string> > & ret );
	bool genotypeEqual(std::string geno,int sample,int snp);
	double getEntropy(std::vector<int> samples, std::vector<int> snps );
	double getInformationInteraction(std::vector<int> samples,std::vector<int> snp);
};

} /* namespace SHEsis */

#endif /* GENEINTERACTION_H_ */
