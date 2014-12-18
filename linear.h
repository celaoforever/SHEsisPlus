/*
 * linear.h
 *
 *  Created on: Dec 18, 2014
 *      Author: ada
 */

#ifndef LINEAR_H_
#define LINEAR_H_
#include <mlpack/core.hpp>
#include <mlpack/methods/linear_regression/linear_regression.hpp>
#include "regression.h"
using namespace mlpack;
using namespace mlpack::regression;
namespace SHEsis {

class linear:public regression {
public:
	linear(std::vector<double>& response, std::vector< std::vector<double> >& _covar,std::vector<double>& snp);
	virtual void regress();
	virtual ~linear();
	virtual void getPvalue();
private:
	arma::vec predictions;
};

} /* namespace SHEsis */

#endif /* LINEAR_H_ */
