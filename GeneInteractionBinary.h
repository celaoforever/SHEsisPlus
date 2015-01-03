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
	std::string snpset;
	double caseEntropy;
	double ctrlEntropy;
	double diff;
};

class GeneInteractionBinary: public GeneInteraction {
public:
	GeneInteractionBinary(boost::shared_ptr<SHEsisData> data);
	virtual ~GeneInteractionBinary();
	virtual void CalGeneInteraction();
	void print();

private:
	std::vector<gxgBinaryRes> res;
	void getNonmissingSample(std::vector<int>& Snp,std::vector<int>& validCase,std::vector<int>& validCtrl);
	void getSingleEntropySum(std::vector<int>& Snp,std::vector<int>& validCase,std::vector<int>& validCtrl,double& CaseSum, double& CtrlSum);
	gxgBinaryRes GetOneSNPCombinationInformationGain(std::vector<int>& Snp);
	void GenerateSNPCombination(int snpnum,std::vector<std::vector<int> >& ret);
	virtual void GetInformationGain(std::vector<int>& Snp,std::vector<std::vector<std::string> >& cp,double& caseGain,double& ctrlGain);
	virtual void GenerateGenotypeCombination(std::vector<int>& Snp,std::vector<std::vector<std::string> > & ret );

};

} /* namespace SHEsis */

#endif /* GENEINTERACTIONBINARY_H_ */
