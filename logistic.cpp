/*
 * logistic.cpp
 *
 *  Created on: Dec 17, 2014
 *      Author: ada
 */

#include "logistic.h"
#include <boost/math/distributions/chi_squared.hpp>
namespace SHEsis {

logistic::logistic(std::vector<double>& response, std::vector< std::vector<double> >& _covar,std::vector<double>& snp):
		regression(response,_covar,snp),
	    optimizerType("lbfgs"),
		positive(2),negative(1),maxIterations(100),tolerance(0.00001)
{};


logistic::~logistic() {
	// TODO Auto-generated destructor stub
}

void logistic::getPvalue(){
	arma::mat ones(this->regressors.size(),1,arma::fill::ones);
	arma::mat X=arma::join_rows(ones,this->regressors);
	arma::mat Xt=X.t();
	arma::mat V(this->responses.size(),this->responses.size(),arma::fill::zeros);
	for (int i=0; i<this->responses.size(); i++)
	{
	    double t = 0;
	    for (int j=0; j<this->regressors.n_cols+1; j++)
	      t += this->coef(j) * X(i,j);
	    double p = 1/(1+exp(-t));
	    V(i,i) = p * (1-p);
	}
	arma::mat S=(Xt*V*X).i();
	for(int i=0;i<p.size();i++){
		this->se(i)=sqrt(S(i,i));
		double Z=this->coef(i)/this->se(i);
		try {
		   boost::math::chi_squared dist(1);
		    this->p(i) = boost::math::cdf(boost::math::complement(dist, Z));
		}catch (...) {
		    this->p(i) = -999;
		  }
	}
}

void logistic::regress(){
	LogisticRegressionFunction lrf(this->regressors,this->responses,this->lambda);
	if(this->optimizerType == "lbfgs"){
		L_BFGS<LogisticRegressionFunction> lbfgsOpt(lrf);
		lbfgsOpt.MaxIterations()=this->maxIterations;
		lbfgsOpt.MinGradientNorm()=this->tolerance;
		LogisticRegression<L_BFGS> lr(lbfgsOpt);
		this->coef=lr.Parameters();
	}else if(this->optimizerType == "sgd"){
		SGD<LogisticRegressionFunction> sgdOpt(lrf);
		sgdOpt.MaxIterations()=this->maxIterations;
		sgdOpt.Tolerance()=this->tolerance;
		sgdOpt.StepSize()=0.01;
		LogisticRegression<SGD> lr(sgdOpt);
		this->coef=lr.Parameters();
	}
	this->getPvalue();
}


} /* namespace SHEsis */
