/*
 * regression.h
 *
 *  Created on: Dec 18, 2014
 *      Author: ada
 */

#ifndef REGRESSION_H_
#define REGRESSION_H_

namespace SHEsis {
#include <mlpack/core.hpp>

using namespace mlpack;

class regression {
public:
	regression(std::vector<double>& response, std::vector< std::vector<double> >& _covar,std::vector<double>& snp);
	virtual ~regression();
	void resetResponse(std::vector<double>& response);
	void resetSnp(std::vector<double>& snp);
	void setLambda(double l){this->lambda=l;};
	virtual void regress()=0;
	arma::vec coef;
	arma::vec p;
	arma::vec se;
protected:
	arma::mat regressors;
	arma::vec responses;

	double lambda;
	virtual void getPvalue()=0;
};

} /* namespace SHEsis */

#endif /* REGRESSION_H_ */
