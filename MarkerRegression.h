/*
 * MarkerRegression.h
 *
 *  Created on: Dec 19, 2014
 *      Author: ionadmin
 */

#ifndef MARKERREGRESSION_H_
#define MARKERREGRESSION_H_
#include "SHEsisData.h"
#include "logistic.h"
#include "linear.h"
#include "CreatHtmlTable.h"
namespace SHEsis {
typedef enum{
	DOMINANT,
	ADDICTIVE,
	RECESSIVE
} DiseaseModel;
struct RegressionRes{
	RegressionRes():p(2),coef(-999),se(-999),allele(-999),nonmissing(-999),HolmP(-999),SidakSSP(-999),
			SidakSDP(-999),BHP(-999),BYP(-999),permutationP(-999),OR(-999),ORlb(-999),ORub(-999){}
	double p;
	double coef;
	double OR;
	double ORlb;
	double ORub;
	double se;
	short allele;
	int nonmissing;
	double HolmP;
	double SidakSSP;
	double SidakSDP;
	double BHP;
	double BYP;
	double permutationP;
	bool operator<(const RegressionRes& res) {
	    if (p < res.p)
	      return true;
	    else
	      return false;
	};
	bool operator>(const RegressionRes& res) {
	    if (p > res.p)
	      return true;
	    else
	      return false;
	};
};
class MarkerRegression {
public:
	std::vector<RegressionRes> vResults;
	boost::shared_ptr<SHEsisData> data;
	MarkerRegression(boost::shared_ptr<SHEsisData> d);
	void setDiseaseModel(DiseaseModel m){this->model=m;};
	void regressAll();
	void setAdjust(bool b){this->adjust=b;};
	virtual ~MarkerRegression();
	void RegressPermutation(std::vector<double>& phenotype);
	std::string reporttxt();
	std::string reporthtml();
	void setPermutation(int p){this->permutation=p;};
private:
	int permutation;
	double codeAllele(int sample,int snp,short allele);
	RegressionRes OneLocusRegression(int snpIdx,short allele);
	int getAlleleCount(int sample, int snp, short allele);
	bool IsGenotypeMissing(int sample, int snp);
	std::vector<short> FindAllele(int snp);
	DiseaseModel model;
	bool adjust;
	boost::shared_ptr<regression> lr;

};
}
#endif /* MARKERREGRESSION_H_ */
