/*
 * GeneInteractionQTL.h
 *
 *  Created on: Jan 10, 2015
 *      Author: ionadmin
 */

#ifndef GENEINTERACTIONQTL_H_
#define GENEINTERACTIONQTL_H_

#include "GeneInteraction.h"
#include "linear.h"
namespace SHEsis {
struct bin{
double start;
double end;
double meanqtl;
double entropy;
std::vector<int> SampleIdx;
};

struct gxgQTLRes{
	//gxgQTLRes():p(-999),nonmissing(-999){}
	gxgQTLRes():p(-999),LowEntropy(-999),HighEntropy(-999),diff(-999),nonmissing(-999){}
	std::string snpset;
	double LowEntropy;
	double HighEntropy;
	double diff;
	double p;
//	double p2;
	double HolmP;
	double SidakSDP;
	double SidakSSP;
	double BHP;
	double BYP;
	int nonmissing;
};

struct qtl2sampleIdx{
	double qtl;
	int idx;
};
class GeneInteractionQTL: public GeneInteraction {
public:
	GeneInteractionQTL(boost::shared_ptr<SHEsisData> data);
	virtual ~GeneInteractionQTL();
	virtual void CalGeneInteraction();
	void print();
private:
	void init();
	gxgQTLRes GetOneSNPCombinationInformationGain(std::vector<int>& Snp);
	void GetNonmissingSamples(std::vector<int>& snp,std::vector<int>& input ,std::vector<int>& output );
	virtual std::string reporttxt();
	virtual std::string reporthtml();
	double getThreshold(std::vector<double> vec,int NumOfBins);

	std::vector<gxgQTLRes> res;
	std::vector<int> lowgroup;
	std::vector<int> highgroup;
};


} /* namespace SHEsis */

#endif /* GENEINTERACTIONQTL_H_ */
