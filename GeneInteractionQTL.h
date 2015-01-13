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
	std::string snpset;
	std::vector<double> entropy; //input
	std::vector<double> qtl; //response
	double coef;
	double stderr;
	double p;
};

class GeneInteractionQTL: public GeneInteraction,linear {
public:
	GeneInteractionQTL(boost::shared_ptr<SHEsisData> data);
	virtual ~GeneInteractionQTL();
	virtual void CalGeneInteraction();
	void setSamplePerBin(int s){this->MinSamplesPerBin=s;};
	void setMaxBin(int s){this->MaxBin=s;};
	void setMinBin(int s){this->MinBin=s;};
private:
	gxgQTLRes GetOneSNPCombinationInformationGain(std::vector<int>& Snp);
	void statIntervalSampleNum(bin& b);
	void GetBins();
	bool ForwardMergeBins(std::list<bin>::iterator& iter);
	bool BackwardMergeBins(std::list<bin>::reverse_iterator& iter);
	std::list<bin> bins;
	std::vector<gxgQTLRes> res;
	int MinSamplesPerBin;
	int MaxBin;
	int MinBin;


};

} /* namespace SHEsis */

#endif /* GENEINTERACTIONQTL_H_ */
