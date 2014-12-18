/*
 * regression.cpp
 *
 *  Created on: Dec 18, 2014
 *      Author: ada
 */

#include "regression.h"

namespace SHEsis {

regression::regression(std::vector<double>& response, std::vector< std::vector<double> >& _covar,std::vector<double>& snp):
		regressors(_covar[0].size(),_covar.size()),lambda(0),
		p(1),coef(1),
		se(1),
		responses(_covar.size())
		{
	// TODO Auto-generated constructor stub
	BOOST_ASSERT(_covar.size() == response.size());
	BOOST_ASSERT(_covar[0].size()!=0);
	for(int i=0;i<_covar.size();i++){
		BOOST_ASSERT(_covar[0].size() == _covar[i].size());
		for(int j=0;j<_covar[0].size();j++){
			regressors(j,i)=_covar[i][j];
		}
		//regressors(i,this->regressors.n_cols-1)=snp[i];
	}
	this->addCovar(snp);

}

void regression::addCovar(std::vector<double> c){
	this->regressors.insert_rows(0,1);
	BOOST_ASSERT(c.size() == this->regressors.n_cols);
	for(int i=0;i<c.size();i++){
		this->regressors(0,i)=c[i];
	}

}

regression::~regression() {
	// TODO Auto-generated destructor stub
}

void regression::resetSnp(std::vector<double>& snp){
	BOOST_ASSERT(snp.size() == this->responses.size());
	for(int i=0;i<this->responses.size();i++){
		this->regressors(i,this->regressors.n_cols-1)=snp[i];
	};
}

void regression::resetResponse(std::vector<double>& response){
	BOOST_ASSERT(this->responses.size() == response.size());
	for(int i=0;i<response.size();i++)
		this->responses[i]=response[i];
}

} /* namespace SHEsis */
