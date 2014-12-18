/*
 * regression.cpp
 *
 *  Created on: Dec 18, 2014
 *      Author: ada
 */

#include "regression.h"

namespace SHEsis {

regression::regression(std::vector<double>& response, std::vector< std::vector<double> >& _covar,std::vector<double>& snp):
		regressors(_covar.size(),_covar[0].size()+1),lambda(0),
		p(this->regressors.n_cols+1),coef(this->regressors.n_cols+1),
		se(this->regressors.n_cols+1),
		responses(_covar.size()) {
	// TODO Auto-generated constructor stub
	BOOST_ASSERT(_covar.size() == response.size());
	BOOST_ASSERT(_covar[0].size()!=0);
	for(int i=0;i<_covar.size();i++){
		if(response[i]==this->positive)
			this->responses(i)=1;
		else if (response[i] == this->negative)
			this->responses(i)=0;
		else
			BOOST_ASSERT(1 == 0);
		this->responses(i)=response[i];
		BOOST_ASSERT(_covar[0].size() == _covar[i].size());
		for(int j=0;j<_covar[0].size();j++){
			regressors(i,j)=_covar[i][j];
		}
		regressors(i,this->regressors.n_cols-1)=snp[i];
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
