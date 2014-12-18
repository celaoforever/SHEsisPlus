/*
 * linear.cpp
 *
 *  Created on: Dec 18, 2014
 *      Author: ada
 */

#include "linear.h"
#include <math.h>
#include <boost/math/distributions/students_t.hpp>
namespace SHEsis {

linear::linear(std::vector<double>& response, std::vector< std::vector<double> >& _covar,std::vector<double>& snp):
		regression(response,_covar,snp){

}
void linear::regress(){
	LinearRegression lr;
	lr.lambda()=this->lambda;
	lr=LinearRegression(this->regressors,this->responses);
	this->coef=lr.Parameters();
	lr.Predict(this->regressors,this->predictions);
	this->getPvalue();
}
void linear::getPvalue(){
	BOOST_ASSERT(this->predictions.size() == this->responses.size());
	double numerator=0;
	for(int i=0;i<this->responses.size();i++){
		numerator+=(this->responses[i]-this->predictions[i])*
				(this->responses[i]-this->predictions[i]);
	}
	numerator=sqrt(numerator/(double)(this->predictions.size()-1));
	arma::vec ave(this->coef.size(),arma::fill::zeros);
	for(int i=1;i<this->coef.size();i++){
		for(int j=0;j<this->responses.size();j++){
			ave(i)+=this->regressors(j,i);
		};
		ave(i)=ave(i)/(double)this->responses.size();
	}
	this->p[0]=-999;
	for(int i=1;i<this->coef.size();i++){//skip the intercept
		double denominator=0;
		int df=0;
		double t=0;
		for(int j=0;j<this->responses.size();j++){
			denominator+=(this->regressors(j,i)-ave(i))*(this->regressors(j,i)-ave(i));
		}
		denominator=sqrt(denominator);
		this->se[i]=numerator/denominator;
		df=this->responses-this->regressors.n_cols;
		t=this->coef[i]/this->se[i];
		try {
		    boost::math::students_t dist(df);
		    this->p[i] = boost::math::cdf(boost::math::complement(dist, t));
		}
		catch (...) {
			this->p[i] = -999;
		}

	}
}

linear::~linear() {
	// TODO Auto-generated destructor stub
}

} /* namespace SHEsis */
