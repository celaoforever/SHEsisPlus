/*
 * logistic_test.cpp
 *
 *  Created on: Dec 17, 2014
 *      Author: ada
 */

#include "logistic.h"
int main(){
	std::vector<double> snp;
	std::vector< std::vector<double> >covar;
	std::vector<double> response;
	//double _snp[6]={2,1,0,1,0,0};
	double _snp[6]={1,1,1,0,0,0};
	double reg[6][3]={
					{39,1,3.2},
					{42,2,3},
					{23,2,4},
					{45,1,3.5},
					{23,1,3},
					{11,2,1.2}
					};
	//double res[6]={1,1,1,2,2,2};
	double res[6]={1,2,2,1,1,2};
	for(int i=0;i<6;i++){
		snp.push_back(_snp[i]);
		response.push_back(res[i]);
	}

	for(int i=0;i<6;i++){
		std::vector<double> t;
		for(int j=0;j<3;j++){
			t.push_back(reg[i][j]);
		}
		covar.push_back(t);
	}

	SHEsis::logistic l;
	l.setCovar(covar);
	l.resetSNP(snp);
	l.setReponse(response);
	l.regress();
	std::cout<<"coef:\n";l.coef.print();
	std::cout<<"se:\n";l.se.print();
	std::cout<<"p:\n";l.p.print();
	return 0;
}
