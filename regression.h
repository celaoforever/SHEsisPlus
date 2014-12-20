/*
 * regression.h
 *
 *  Created on: Dec 18, 2014
 *      Author: ada
 */

#ifndef REGRESSION_H_
#define REGRESSION_H_
#include <mlpack/core.hpp>
using namespace mlpack;
namespace SHEsis {


class regression {
public:
	regression();
	virtual ~regression();
	void resetResponse(std::vector<double>& response);
	void modifyCovar(std::vector<double>& snp,int row);
	void setLambda(double l){this->lambda=l;};
	void pushCovar(std::vector<double>& c);
	void resetSNP(std::vector<double>& c);
	virtual void regress(){};
	void setCovar(std::vector< std::vector<double> >& _covar);
	virtual void setReponse(std::vector<double>& response){};
	arma::vec coef;
	arma::vec p;
	arma::vec se;
protected:
	arma::mat regressors;
	arma::vec responses;
	bool SNPAdded;
	double lambda;
	virtual void getPvalue(){};
};

} /* namespace SHEsis */

#endif /* REGRESSION_H_ */
