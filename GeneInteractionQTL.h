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
	gxgQTLRes():p(-999),nonmissing(-999),numBin(-999){}
	std::string snpset;
	double LowEntropy;
	double HighEntropy;
	double diff;
	double p;
	double p2;
	double HolmP;
	double SidakSDP;
	double SidakSSP;
	double BHP;
	double BYP;
	int nonmissing;
	//added
std::vector<double> entropy; //input
std::vector<double> qtl; //response
double rankcorrel;
double rankporigin;
double rankp;
double rankp2;
int numBin;
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
	//added
	virtual void setBinNum(int s){this->NumBin=s;};
	virtual void setMinBin(int s){this->MinBin=s;};
	virtual void setSamplePerBin(int s){this->MinSamplesPerBin=s;};
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

	//added
	std::list<bin> bins;
	int MinSamplesPerBin;
	int NumBin;
	int MinBin;
	boost::shared_ptr<linear> lr;
	double getIntervalQTL(bin& b);
	void randomshuffle(std::vector<qtl2sampleIdx>& samples);
	int UpdateBins2(std::vector<int>& snp,std::vector<qtl2sampleIdx>& validSamples);
	void GetNonmissingSamples2(std::vector<int>& snp,std::vector<qtl2sampleIdx>& output );
	gxgQTLRes GetOneSNPCombinationInformationGain2(std::vector<int>& Snp);
	void GetRankCorrelation(std::vector<double> entropy, double& correl,double& p);
	void MergeBins(std::vector<int>& Snp,std::list<bin>& bins,int count,std::vector<double>& entropy);
	bool ForwardMergeBins(std::list<bin>::iterator& iter,std::vector<int>& Snps);
	std::list<bin>::iterator findMin(std::list<bin>& bins,double& mindiff);

};


} /* namespace SHEsis */

#endif /* GENEINTERACTIONQTL_H_ */
