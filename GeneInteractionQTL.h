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
	std::vector<int> SampleIdx;
};

struct gxgQTLRes{
	gxgQTLRes():coef(-999),se(-999),p(-999),nonmissing(-999),numBin(-999),correl(-999){}
	std::string snpset;
	std::vector<double> entropy; //input
	std::vector<double> qtl; //response
	double diff;
	double correl;
	double coef;
	double rankcorrel;
	double rankporigin;
	double rankp;
	double rankp2;
	double se;
	double p;
	double p2;
	double HolmP;
	double SidakSDP;
	double SidakSSP;
	double BHP;
	double BYP;
	int numBin;
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
	virtual void setBinNum(int s){this->NumBin=s;};
	virtual void  setMinBin(int s){this->MinBin=s;};
	virtual void setSamplePerBin(int s){this->MinSamplesPerBin=s;};
	void seta(double a){this->a=a;}
	void print();
private:
	double getP(std::vector<double> permutated, double origin);
	gxgQTLRes GetOneSNPCombinationInformationGain(std::vector<int>& Snp);
	double getIntervalQTL(bin& b);
	void statIntervalSampleNum(bin& b);
	void statIntervalSampleNum(bin& b,std::vector<int>& samples);
	int  UpdateBins(std::vector<int>& snp);
	int  UpdateBins2(std::vector<int>& snp,std::vector<qtl2sampleIdx>& validsample);
	void randomshuffle(std::vector<qtl2sampleIdx>& validsample);
	double GetCorrelationCoeff(std::vector<double> entropy,std::vector<double> qtl);
	void GetRankCorrelation(std::vector<double> entropy, double& correl,double& p);
	void GetBins();
	void GetNonmissingSamples(std::vector<int>& snp,std::vector<int>& output );
	void GetNonmissingSamples2(std::vector<int>& snp,std::vector<qtl2sampleIdx>& output );
	bool ForwardMergeBins(std::list<bin>::iterator& iter);
	bool BackwardMergeBins(std::list<bin>::reverse_iterator& iter);	std::list<bin> bins;
	virtual std::string reporttxt();
	virtual std::string reporthtml();
	std::vector<gxgQTLRes> res;
	boost::shared_ptr<linear> lr;
	int MinSamplesPerBin;
	int NumBin;
	int MinBin;
	double a;
};


} /* namespace SHEsis */

#endif /* GENEINTERACTIONQTL_H_ */
