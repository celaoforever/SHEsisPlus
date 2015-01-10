/*
 * GeneInteractionBinary.h
 *
 *  Created on: Jan 1, 2015
 *      Author: ionadmin
 */

#ifndef GENEINTERACTIONBINARY_H_
#define GENEINTERACTIONBINARY_H_
#define GXG_CASE 2
#define GXG_CTRL 1
#include "GeneInteraction.h"
namespace SHEsis {


struct gxgBinaryRes{
	gxgBinaryRes():permutatedDiffMean(0),permutatedDiffVar(0),p(-999){};
	std::string snpset;
	double caseEntropy;
	double ctrlEntropy;
	double diff;
	double caseLambda;
	double ctrlLambda;
	std::vector<double> permutatedCaseEntropy;
	std::vector<double> permutatedCtrlEntropy;
	double permutatedDiffMean;
	double permutatedDiffVar;
	double p;
};

class GeneInteractionBinary: public GeneInteraction {
public:
	GeneInteractionBinary(boost::shared_ptr<SHEsisData> data);
	virtual ~GeneInteractionBinary();
	virtual void CalGeneInteraction();
	void print();
	void setPermutation(int p){this->permutation=p;};
private:
	int permutation;
	std::vector<gxgBinaryRes> res;
	gxgBinaryRes GetOneSNPCombinationInformationGain2(std::vector<int>& Snp);
	void getNonmissingSample(std::vector<int>& Snp,std::vector<int>& validCase,std::vector<int>& validCtrl);
	//below not used in new version
	void getSingleEntropySum(std::vector<int>& Snp,std::vector<int>& validCase,std::vector<int>& validCtrl,double& CaseSum, double& CtrlSum);
	gxgBinaryRes GetOneSNPCombinationInformationGain(std::vector<int>& Snp);
	void GetInformationGain(std::vector<int>& Snp,std::vector<std::vector<std::string> >& cp,double& caseGain,double& ctrlGain);
};

} /* namespace SHEsis */

#endif /* GENEINTERACTIONBINARY_H_ */
