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
	std::vector<int> SampleIdx;
};

struct gxgQTLRes{
	gxgQTLRes():coef(-999),se(-999),p(-999),nonmissing(-999),numBin(-999){}
	std::string snpset;
	std::vector<double> entropy; //input
	std::vector<double> qtl; //response
	double coef;
	double se;
	double p;
	int numBin;
	int nonmissing;
};

class GeneInteractionQTL: public GeneInteraction {
public:
	GeneInteractionQTL(boost::shared_ptr<SHEsisData> data);
	virtual ~GeneInteractionQTL();
	virtual void CalGeneInteraction();
	void setSamplePerBin(int s){this->MinSamplesPerBin=s;};
	void setMaxBin(int s){this->MaxBin=s;};
	void setMinBin(int s){this->MinBin=s;};
	void print();
private:
	gxgQTLRes GetOneSNPCombinationInformationGain(std::vector<int>& Snp);
	double getIntervalQTL(bin& b);
	void statIntervalSampleNum(bin& b);
	void statIntervalSampleNum(bin& b,std::vector<int>& samples);
	int  UpdateBins(std::vector<int>& snp);
	void GetBins();
	void GetNonmissingSamples(std::vector<int>& snp,std::vector<int>& output );
	bool ForwardMergeBins(std::list<bin>::iterator& iter);
	bool BackwardMergeBins(std::list<bin>::reverse_iterator& iter);
	std::list<bin> bins;
	std::vector<gxgQTLRes> res;
	boost::shared_ptr<linear> lr;
	int MinSamplesPerBin;
	int MaxBin;
	int MinBin;
};

} /* namespace SHEsis */

#endif /* GENEINTERACTIONQTL_H_ */
