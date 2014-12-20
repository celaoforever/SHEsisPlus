/*
 * MarkerRegression_test.cpp
 *
 *  Created on: Dec 20, 2014
 *      Author: ionadmin
 */
#include "MarkerRegression.h"
using namespace SHEsis;
using namespace std;

int main(){
//	double snp[6][2][2]={
//			{{1,1},{1,2},{2,2},{1,2},{2,2},{2,2}},
//			{{1,2},{1,2},{1,2},{1,1},{1,1},{1,1}}
//	};
//	double snp[6][2][2]={
//			{{1,1},{1,2}},
//			{{1,2},{1,2}},
//			{{2,2},{1,2}},
//			{{1,2},{1,1}},
//			{{2,2},{1,1}},
//			{{2,2,},{1,1}}
//
//	};
	string snp[6][2][2]={
			{{"1","1"},{"1","2"}},
			{{"1","2"},{"1","2"}},
			{{"2","2"},{"1","2"}},
			{{"1","2"},{"1","1"}},
			{{"2","2"},{"1","1"}},
			{{"2","2"},{"1","1"}}

	};
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
	int sampleNum = 6;
	int snpNum = 2;
	int ploidy = 2;
	boost::shared_ptr<SHEsisData> data(new SHEsisData(sampleNum,snpNum,ploidy));
	//covar
	for(int i=0;i<sampleNum;i++){
		std::vector<double> row;
		for(int j=0;j<3;j++){
			row.push_back(reg[i][j]);
		}
		data->covar.push_back(row);
	}
	//genotype and phenotype
	for(int i=0;i<sampleNum;i++){
		data->vLabel.push_back((SampleStatus)((int)res[i]));
		for(int j=0;j<snpNum;j++){
			for(int p=0;p<ploidy;p++){
				data->mGenotype[i][j][p]=data->GetAlleleCode(snp[i][j][p]);
			}
		}
	}
	cout<<"logistic:\n";
	MarkerRegression mr(data);
	mr.setAdjust(true);
	mr.regressAll();
	//cout<<mr.reporthtml();
	cout<<mr.reporttxt();
	data->vLabel.clear();
	double res2[6]={4,2,0,2,0,0};
	for(int i=0;i<sampleNum;i++){
		data->vQuantitativeTrait.push_back(res2[i]);
	}
	cout<<"linear:\n";
	mr.regressAll();
	cout<<mr.reporttxt();


	return 0;
}



