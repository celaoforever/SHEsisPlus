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

linear::linear():
		regression(){};

void linear::setReponse(std::vector<double>& response){
	this->responses.resize(response.size());
	for(int i=0;i<response.size();i++){
		this->responses(i)=response[i];
	}
}
void linear::regress(){
//	std::cout<<"regressors:\n";
//	regressors.print();
//	std::cout<<"responses:\n";
//	responses.print();
	BOOST_ASSERT(this->regressors.n_cols>0 && this->regressors.n_rows>0 && this->responses.size()>0);
	BOOST_ASSERT(this->regressors.n_cols==this->responses.size());
	LinearRegression lr;
	lr.Lambda()=this->lambda;
	lr=LinearRegression(this->regressors,this->responses);
	this->coef=lr.Parameters();
	lr.Predict(this->regressors,this->predictions);
	this->getPvalue();
//	std::cout<<"coef:\n";coef.print();
}
void linear::getPvalue(){
	BOOST_ASSERT(this->predictions.size() == this->responses.size());
	this->p.resize(this->regressors.n_rows+1);
	this->se.resize(this->regressors.n_rows+1);
	double numerator=0;
	for(int i=0;i<this->responses.size();i++){
		numerator+=(this->responses(i)-this->predictions(i))*
				(this->responses(i)-this->predictions(i));
	}
	numerator=sqrt(numerator/(double)(this->predictions.size()-this->regressors.n_rows));
	arma::vec ave(this->coef.size(),arma::fill::zeros);
	for(int i=1;i<this->coef.size();i++){
		for(int j=0;j<this->responses.size();j++){
			ave(i)+=this->regressors(i-1,j);
		};
		ave(i)=ave(i)/(double)this->responses.size();
	}
	this->p(0) = 1;
	for(int i=1;i<this->coef.size();i++){
		double denominator=0;
		int df=0;
		double t=0;
		for(int j=0;j<this->responses.size();j++){
			denominator+=(this->regressors(i-1,j)-ave(i))*(this->regressors(i-1,j)-ave(i));
		}
		denominator=sqrt(denominator);
		this->se[i]=numerator/denominator;
		df=this->responses.size()-this->regressors.n_rows;
		t=this->coef[i]/this->se[i];
		try {
			t=t>0?t:(-1*t);
		    boost::math::students_t dist(df);
//		    if(i == 1)
//		    	std::cout<<"t="<<t<<"\n";
		    this->p(i) = 2*boost::math::cdf(boost::math::complement(dist, t));
		}
		catch (...) {
			this->p(i) = -999;
		}

	}
}

linear::~linear() {
	// TODO Auto-generated destructor stub
}

} /* namespace SHEsis */
