/*
 * regression.cpp
 *
 *  Created on: Dec 18, 2014
 *      Author: ada
 */

#include "regression.h"

namespace SHEsis {

regression::regression():
		lambda(0),SNPAdded(false)
{
	// TODO Auto-generated constructor stub

	//this->pushCovar(snp);
}
void regression::setCovar(std::vector< std::vector<double> >& _covar){
	BOOST_ASSERT(_covar[0].size()!=0);
	this->regressors.resize(_covar[0].size(),_covar.size());
	//calculate average
	std::vector<double> ave(_covar[0].size(),0);
	std::vector<double> count(_covar[0].size(),0);
	for(int i=0;i<_covar.size();i++){
		BOOST_ASSERT(_covar[0].size() == _covar[i].size());
		for(int j=0;j<_covar[0].size();j++){
			if(_covar[i][j] != DBL_MAX){
				count[j]++;
				ave[j]+=_covar[i][j];
			}
		}
	}
	for(int j=0;j<_covar[0].size();j++){
		ave[j]/=(double)count[j];
	}
	//fill missing data with average
	for(int i=0;i<_covar.size();i++){
		BOOST_ASSERT(_covar[0].size() == _covar[i].size());
		for(int j=0;j<_covar[0].size();j++){
			if(_covar[i][j] == DBL_MAX){
				regressors(j,i)=ave[j];
			}else{
				regressors(j,i)=_covar[i][j];
			}
		}
	}
}

void regression::resetSNP(std::vector<double>& c){
	if(this->SNPAdded){
		this->modifyCovar(c,0);
	}else{
		this->SNPAdded=true;
		this->pushCovar(c);
	}
}

void regression::pushCovar(std::vector<double>& c){
	if(this->regressors.n_cols!= 0){
		BOOST_ASSERT(c.size() == this->regressors.n_cols);
		this->regressors.insert_rows(0,1);
	}else{
		this->regressors.resize(1,c.size());
	}
	for(int i=0;i<c.size();i++){
		this->regressors(0,i)=c[i];
	}
}

regression::~regression() {
	// TODO Auto-generated destructor stub
}

void regression::modifyCovar(std::vector<double>& snp,int row){
	BOOST_ASSERT(snp.size() == this->responses.size());
	for(int i=0;i<this->responses.size();i++){
		this->regressors(row,i)=snp[i];
	};
}

void regression::resetResponse(std::vector<double>& response){
	this->responses.resize(response.size());
	for(int i=0;i<response.size();i++)
		this->responses[i]=response[i];
}

} /* namespace SHEsis */
