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
	for(int i=0;i<response.size();i++){
		this->responses(i)=response[i];
	}
}
void linear::regress(){
	LinearRegression lr;
	lr.Lambda()=this->lambda;
	lr=LinearRegression(this->regressors,this->responses);
	this->coef=lr.Parameters();
	std::cout<<"coeff:\n";this->coef.print();
	lr.Predict(this->regressors,this->predictions);
	this->getPvalue();
}
void linear::getPvalue(){
	BOOST_ASSERT(this->predictions.size() == this->responses.size());
	this->p.resize(this->regressors.n_rows);
	this->se.resize(this->regressors.n_rows);
	double numerator=0;
	for(int i=0;i<this->responses.size();i++){
		numerator+=(this->responses(i)-this->predictions(i))*
				(this->responses(i)-this->predictions(i));
	}
	std::cout<<"prediction:\n";this->predictions.print();
	std::cout<<"response:\n";this->responses.print();
	numerator=sqrt(numerator/(double)(this->predictions.size()-2));
	arma::vec ave(this->coef.size()-1,arma::fill::zeros);
	for(int i=0;i<this->coef.size()-1;i++){
		for(int j=0;j<this->responses.size();j++){
			ave(i)+=this->regressors(i,j);
		};
		ave(i)=ave(i)/(double)this->responses.size();
	}
	for(int i=0;i<this->coef.size()-1;i++){//skip the intercept
		double denominator=0;
		int df=0;
		double t=0;
		for(int j=0;j<this->responses.size();j++){
			denominator+=(this->regressors(i,j)-ave(i))*(this->regressors(i,j)-ave(i));
		}
		denominator=sqrt(denominator);
		this->se[i]=numerator/denominator;
		df=this->responses.size()-this->regressors.n_rows;
		t=this->coef[i+1]/this->se[i];
		std::cout<<t<<",";
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
