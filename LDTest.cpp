/*
 * LDTest.cpp
 *
 *  Created on: Sep 17, 2014
 *      Author: ada
 */

#include "LDTest.h"
#include "HaplotypeDiploid.h"
#include "Haplotype.h"
namespace SHEsis {

LDTest::LDTest(boost::shared_ptr<SHEsisData> data):
		data(data),ldtype(LD_IN_BOTH),
		res(boost::extents[this->data->getSnpNum()][this->data->getSnpNum()]) {
	// TODO Auto-generated constructor stub

}

LDTest::~LDTest() {
	// TODO Auto-generated destructor stub
}
//type=0: case, type=1:ctrl, type=2 both
double LDTest::TwoLociLDTest(int snp1,int snp2,LD_TYPE type){
	double normalizedD=0;
	std::vector<short> mask;
	for(int i=0;i<this->data->getSnpNum();i++){
		if(i == snp1 || i == snp2)
			mask.push_back(1);
		else
			mask.push_back(0);
	}
	BOOST_ASSERT(mask.size()==this->data->getSnpNum());
	if(this->data->getNumOfChrSet()<= 2){
		this->hp.reset(new HaplotypeDiploid(this->data,2,mask));
	}else{
		this->hp.reset(new Haplotype(this->data,2,mask));
	}
	hp->startHaplotypeAnalysis();

	for(int i=0;i<hp->Results.haplotypes.size();i++){
		double hapfreq=0;
		double allelefreq1=0;
		double allelefreq2=0;
		int allele=0;
		double d=0;
		double dmax=0;
		switch(type){
		case LD_IN_CASE:
			hapfreq=(double)hp->Results.CaseCount[i]/(double)(2*this->data->getCaseNum());
			allele=hp->Results.haplotypes[i][0];
			allelefreq1=(double)data->vLocusInfo[snp1].CaseAlleleCount[allele]/
					(double)(this->data->getNumOfChrSet()*this->data->getCaseNum());
			allele=hp->Results.haplotypes[i][1];
			allelefreq2=(double)data->vLocusInfo[snp2].CaseAlleleCount[allele]/
					(double)(this->data->getNumOfChrSet()*this->data->getCaseNum());
			break;
		case LD_IN_CTRL:
			hapfreq=(double)hp->Results.ControlCount[i]/(double)(this->data->getNumOfChrSet()*this->data->getControlNum());
			allele=hp->Results.haplotypes[i][0];
			allelefreq1=(double)data->vLocusInfo[snp1].ControlAlleleCount[allele]/
					(double)(this->data->getNumOfChrSet()*this->data->getControlNum());
			allele=hp->Results.haplotypes[i][1];
			allelefreq2=(double)data->vLocusInfo[snp2].ControlAlleleCount[allele]/
					(double)(this->data->getNumOfChrSet()*this->data->getCaseNum());
			break;
		case LD_IN_BOTH:
			hapfreq=(double)(hp->Results.CaseCount[i]+hp->Results.ControlCount[i])/(double)(this->data->getNumOfChrSet()*this->data->getSampleNum());
			allele=hp->Results.haplotypes[i][0];
			allelefreq1=(double)data->vLocusInfo[snp1].BothAlleleCount[allele]/
					(double)(this->data->getNumOfChrSet()*this->data->getSampleNum());
			allele=hp->Results.haplotypes[i][1];
			allelefreq2=(double)data->vLocusInfo[snp2].BothAlleleCount[allele]/
					(double)(this->data->getNumOfChrSet()*this->data->getSampleNum());
			break;
		default:
			std::cout<<"***ERROR:no such ld analysis type\n";
			exit(-1);
		};
		d=hapfreq-allelefreq1*allelefreq2;
		double pq=allelefreq1*allelefreq2;
		double p1q1=(1-allelefreq1)*(1-allelefreq2);
		double p1q=allelefreq1*(1-allelefreq2);
		double q1p=allelefreq2*(1-allelefreq1);
		if(d<0){
			dmax=pq<p1q1?pq:p1q1;
		}else{
			dmax=p1q<q1p?p1q:q1p;
		}
		normalizedD+=allelefreq1*allelefreq2*abs(d/dmax);
	}
	return normalizedD;
}


void LDTest::AllLociLDtest(){
	for(int i=0;i<this->data->getSnpNum();i++){
		for(int j=i+1;j<this->data->getSnpNum();j++){
			this->res[i][j]=this->TwoLociLDTest(i,j,this->ldtype);
		}
	}
}
void LDTest::DrawLDMap(){
	double snpnum=(double)this->data->getSnpNum();
	double height,width;
	int recnum;
	recnum=snpnum-1;
	if(snpnum<20)
		height=500;
	else
		height=snpnum*25;
	width=(height-230)/snpnum*2*(snpnum-1)+100;
	double sidelength=(height-230)/snpnum*2;
	double insidelength=sidelength-2;
	double from_x=(double)width/2-(double)recnum*(double)sidelength/2;
	this->ldmap=BMP_new(width,height);
	BMP_clear(this->ldmap, RGB_GRAY);

	for(int j=0;j<recnum;j++){
		int count=0;
		for(int i=1;i<recnum-j;i++){
			int step=i*sidelength;
			int begin_x=from_x+sidelength/2*j;
			int begin_y=200+sidelength/2*j;
			double score=rand()%50+50;
			int x1=round(begin_x+step-sidelength/2);
			int y1=round(begin_y+1);
			int x2=round(begin_x+step-1);
			int y2=round(begin_y+sidelength/2);
			int x3=round(begin_x+step-sidelength/2);
			int y3=round(begin_y+1+sidelength);
			int x4=round(begin_x+step-sidelength+1);
			int y4=round(begin_y+sidelength/2);
			BMP_line (this->ldmap, x1,y1,x2,y2, RGB_CYAN);
			BMP_line (this->ldmap, x2,y2,x3,y3, RGB_CYAN);
			BMP_line (this->ldmap, x3,y3,x4,y4, RGB_CYAN);
			BMP_line (this->ldmap, x4,y4,x1,y1, RGB_CYAN);
			std::stringstream ss;
			ss<<score;
			int strx=begin_x+step-sidelength/2-7;
			int stry=begin_y+sidelength/2-7;
			BMP_draw_string(this->ldmap,ss.str().c_str(),strx,stry,RGB_PURPLE);
		}
	}
	BMP_write (this->ldmap, "ldmap.bmp");


}


} /* namespace SHEsis */
