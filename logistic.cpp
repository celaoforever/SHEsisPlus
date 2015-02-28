/*
 * logistic.cpp
 *
 *  Created on: Dec 17, 2014
 *      Author: ada
 */

#include "logistic.h"
#include <boost/math/distributions/chi_squared.hpp>
namespace SHEsis {

logistic::logistic():
		regression(),
	    optimizerType("lbfgs"),
		positive(2),negative(1),maxIterations(0),tolerance(0.00000001)
{


};


logistic::~logistic() {
	// TODO Auto-generated destructor stub
}
void logistic::setReponse(std::vector<double>& response){
	this->responses.resize(response.size());
	for(int i=0;i<response.size();i++){
		if(response[i]==this->positive)
			this->responses(i)=1;
		else if (response[i] == this->negative)
			this->responses(i)=0;
		else
			BOOST_ASSERT(1 == 0);
	}

}

void logistic::getPvalue(){
	this->p.resize(this->regressors.n_rows+1);
	this->se.resize(this->regressors.n_rows+1);
	arma::mat ones(1,this->regressors.n_cols,arma::fill::ones);
	arma::mat X=arma::join_cols(ones,this->regressors);
	arma::mat Xt=X.t();
	arma::mat V(this->responses.size(),this->responses.size(),arma::fill::zeros);
	for (int i=0; i<this->responses.size(); i++)
	{
	    double t = 0;
	    for (int j=0; j<this->coef.size(); j++)
	      t += this->coef(j) * Xt(i,j);
	    double p = 1/(1+exp(-t));
	    V(i,i) = p * (1-p);
	}
	arma::mat S;
	try{
		S=(X*V*Xt).i();
	}catch(...){
		for(int i=0;i<p.size();i++){
			this->p(i) = -999;
			this->se(i) = -999;
			this->coef(i) = -999;
		}
		return;
	}
	for(int i=0;i<p.size();i++){
		this->se(i)=sqrt(S(i,i));
		double Z=this->coef(i)/this->se(i);
		try {
		   boost::math::chi_squared dist(1);
		    this->p(i) =2* boost::math::cdf(boost::math::complement(dist, Z*Z));
		}catch (...) {
		    this->p(i) = -999;
		  }
	}
}

void logistic::regress(){
//	std::cout<<"regressors:\n";
//	regressors.print();
//	std::cout<<"responses:\n";
//	responses.print();
	BOOST_ASSERT(this->regressors.n_cols>0 && this->regressors.n_rows>0 && this->responses.size()>0);
	BOOST_ASSERT(this->regressors.n_cols==this->responses.size());
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
//	std::cout<<"coef:\n";
//	coef.print();

}


} /* namespace SHEsis */
