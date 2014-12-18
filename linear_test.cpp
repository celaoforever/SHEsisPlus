/*
 * linear_test.cpp
 *
 *  Created on: Dec 18, 2014
 *      Author: ada
 */

#include "linear.h"
/*
2,39,1,3.2,4
1,42,2,3,1
0,23,2,4,1
1,45,1,3.5,2
0,23,1,3,2
0,11,2,1.2,2
 */
int main(){
	std::vector<double> snp;
	std::vector< std::vector<double> >covar;
	std::vector<double> response;
	double _snp[6]={2,1,0,1,0,0};
	double reg[6][3]={
					{39,1,3.2},
					{42,2,3},
					{23,2,4},
					{45,1,3.5},
					{23,1,3},
					{11,2,1.2}
					};
	double res[6]={4,2,0,2,0,0};
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

	SHEsis::linear l(response,covar,snp);
	l.regress();
	std::cout<<"coef:";
	for(int i=0;i<l.coef.size();i++){
		std::cout<<l.coef(i)<<",";
	}
	std::cout<<"\nse:";
	for(int i=0;i<l.se.size();i++){
		std::cout<<l.se(i)<<",";
	}
	std::cout<<"\np:";
	for(int i=0;i<l.p.size();i++){
		std::cout<<l.p(i)<<",";
	}
	return 0;
}
