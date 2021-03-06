/*
 * logistic.h
 *
 *  Created on: Dec 17, 2014
 *      Author: ada
 */

#ifndef LOGISTIC_H_
#define LOGISTIC_H_
#include <mlpack/methods/logistic_regression/logistic_regression.hpp>
#include <mlpack/core/optimizers/sgd/sgd.hpp>
#include "regression.h"
using namespace mlpack;
using namespace mlpack::regression;
using namespace mlpack::optimization;
namespace SHEsis {


struct logisticRes{

};

class logistic:public regression {
public:
	logistic();
	virtual ~logistic();
	virtual void regress();
	void setReponse(std::vector<double>& response);
	void setClass(double pos,double neg){this->positive=pos;this->negative=neg;}
	void setTolerance(double t){this->tolerance=t;}
	void setMaxInterations(int v){this->maxIterations=v;}
	void setOptimizerType(std::string& t){
		BOOST_ASSERT(t=="sgd" || t == "lbfgs");
		this->optimizerType=t;}

protected:
	virtual void getPvalue();
	std::string optimizerType;
	int positive;
	int negative;
	int maxIterations;
	double tolerance;
};

} /* namespace SHEsis */

#endif /* LOGISTIC_H_ */
