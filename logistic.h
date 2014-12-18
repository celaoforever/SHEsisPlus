/*
 * logistic.h
 *
 *  Created on: Dec 17, 2014
 *      Author: ada
 */

#ifndef LOGISTIC_H_
#define LOGISTIC_H_
#include "SHEsisData.h"
#include <mlpack/core.hpp>
#include <mlpack/methods/logistic_regression/logistic_regression.hpp>
#include <mlpack/core/optimizers/sgd/sgd.hpp>

using namespace std;
using namespace mlpack;
using namespace mlpack::regression;
using namespace mlpack::optimization;
namespace SHEsis {
//typedef enum{
//	DOMINANT,
//	ADDICTIVE,
//	RECESSIVE
//} Disease_model;

struct logisticRes{

};

class logistic {
public:
	logistic(std::vector<double>& response, std::vector< std::vector<double> >& _covar,std::vector<double>& snp);
	virtual ~logistic();
//	void setDiseaseModel(Disease_model m){this->model=m;}
	void resetResponse(std::vector<double>& response);
	void resetSnp(std::vector<double>& snp);
	void logisticRegression();//return the cofficient
	void setClass(double pos,double neg){this->positive=pos;this->negative=neg;}
	void setLambda(double l){this->lambda=l;};
	void setTolerance(double t){this->tolerance=t;}
	void setMaxInterations(int v){this->maxIterations=v;}
	void setOptimizerType(std::string& t){
//		std::string t1="sgd";
//		std::string t2="lbfgs";
		BOOST_ASSERT(t=="sgd" || t == "lbfgs");
//		BOOST_ASSERT(std::strcmp(t.c_str(),t1.c_str()) || std::strcmp(t.c_str(),t2.c_str()));
		this->optimizerType=t;}

private:
	arma::mat regressors;
	arma::vec responses;
	arma::vec coef;
	arma::vec p;
	arma::vec se;
	std::string optimizerType;
	double lambda;
	int positive;
	int negative;
	int maxIterations;
	double tolerance;
	void getPvalue();


};

} /* namespace SHEsis */

#endif /* LOGISTIC_H_ */
