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
	gxgBinaryRes():permutatedDiffMean(0),permutatedDiffVar(0),p(-999),nonmissing(-999){};
	std::string snpset;
	double caseEntropy;
	double ctrlEntropy;
	double diff;
	int nonmissing;
	std::vector<double> permutatedCaseEntropy;
	std::vector<double> permutatedCtrlEntropy;
	double permutatedDiffMean;
	double permutatedDiffVar;
	double p;
//	double p2;
	double HolmP;
	double SidakSDP;
	double SidakSSP;
	double BHP;
	double BYP;
};

class GeneInteractionBinary: public GeneInteraction {
public:
	GeneInteractionBinary(boost::shared_ptr<SHEsisData> data);
	virtual ~GeneInteractionBinary();
	virtual void CalGeneInteraction();
	void print();
private:
	virtual std::string reporttxt();
	virtual std::string reporthtml();
	std::vector<gxgBinaryRes> res;
	gxgBinaryRes GetOneSNPCombinationInformationGain(std::vector<int>& Snp);
	void getNonmissingSample(std::vector<int>& Snp,std::vector<int>& validCase,std::vector<int>& validCtrl);
};

} /* namespace SHEsis */

#endif /* GENEINTERACTIONBINARY_H_ */
