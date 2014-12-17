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
#include "logistic_regression.hpp"
#include <mlpack/core/optimizers/sgd/sgd.hpp>

using namespace std;
using namespace mlpack;
using namespace mlpack::regression;
using namespace mlpack::optimization;
namespace SHEsis {
typedef enum{
	DOMINANT,
	ADDICTIVE,
	RECESSIVE
} Disease_model;

struct logisticRes{

};

class logistic {
public:
	boost::shared_ptr<SHEsisData> data;
	logistic(boost::shared_ptr<SHEsisData> data, std::vector< std::vector<double> > _covar):data(data),model(ADDICTIVE),
			covar(_covar.size(),_covar[0].size()),
			covar(_covar.size()){
		BOOST_ASSERT(_covar.size() == this->data->getSampleNum());
		BOOST_ASSERT(_covar[0].size()!=0);
		for(int i=0;i<_covar.size();i++){
			BOOST_ASSERT(_covar[0].size() == _covar[i].size());
			for(int j=0;j<_covar[0].size();j++){
				covar(i,j)=_covar[i][j];
			}
		}
	};
	virtual ~logistic();
	void setDiseaseModel(Disease_model m){this->model=m;}
	logisticRes singleSnpLogisticRegression(int snpIdx,short testAllele);
	arma::mat& logisticRegression(arma::mat& reg,arma::mat &res);//return the cofficient
	double getPvalue(arma::mat& reg,arma::mat& para);


private:
	arma::mat regressors;
	arma::mat covar;
	arma::mat responses;
	Disease_model model;



};

} /* namespace SHEsis */

#endif /* LOGISTIC_H_ */
