/*
 * HaplotypeBase.cpp
 *
 *  Created on: Sep 20, 2014
 *      Author: ionadmin
 */
#include "HaplotypeBase.h"
#include "CreatHtmlTable.h"
template std::string convert2string<double>(double v);
namespace SHEsis{

void HaplotypeBase::AssociationTest(){

	int haploNum=this->Results.haplotypes.size();
	this->Results.singleHap.resize(haploNum);
    int nrow=2;
	double expect = -1.0;
	double percnt = 100.0;
	double emin = 0;
	double pre = 0, prt = 0;
	int ws = 300000;

	int validHap=0;
    double* contigency= new double[4];
	//single hapotype statistics
	for(int i=0;i<haploNum;i++){
		double freq=((double)this->Results.ControlCount[i]+(double)this->Results.CaseCount[i])/((double)this->data->getSampleNum()*(double)this->data->getNumOfChrSet());
		if(freq<this->freqthreshold)
			continue;
		validHap++;
		contigency[0]=this->Results.ControlCount[i];
		contigency[1]=this->Results.CaseCount[i];
		contigency[2]=this->data->getControlNum()*this->data->getNumOfChrSet()-contigency[0];
		contigency[3]=this->data->getCaseNum()*this->data->getNumOfChrSet()-contigency[1];
		  try{
			  fexact(&nrow, &nrow, contigency, &nrow, &expect, &percnt, &emin, &prt, &pre, &ws);
			  this->Results.singleHap[i].fisherp=pre;
		  }catch(std::runtime_error &){
			  this->Results.singleHap[i].fisherp=-1;
		  }
		  PearsonChiSquareTest(contigency,nrow,nrow,this->Results.singleHap[i].chisquare,this->Results.singleHap[i].pearsonp);
		  if(( 0 != contigency[3] && 0!= contigency[0]
		      && 0 !=contigency[1] && 0 != contigency[2])){
			  this->Results.singleHap[i].OR=(contigency[1]*contigency[2]/(contigency[0]*contigency[3]));
			  double v=1/contigency[1]+1/contigency[2]+1/contigency[3]+1/contigency[0];
			  this->Results.singleHap[i].orlow=this->Results.singleHap[i].OR*exp(-1.96*sqrt(v));
			  this->Results.singleHap[i].orUp=this->Results.singleHap[i].OR*exp(1.96*sqrt(v));
		  }else
		  {
			  this->Results.singleHap[i].OR=-1;
			  this->Results.singleHap[i].orlow=-1;
			  this->Results.singleHap[i].orUp=-1;
		  }


	}

	delete[] contigency;
//    std::cout<<"contigency:\n";
	if(validHap==0)
		return;
	contigency= new double[2*validHap];
    int idx=0;
    for(int i=0;i<haploNum;i++){
		double freq=((double)this->Results.ControlCount[i]+(double)this->Results.CaseCount[i])/((double)this->data->getSampleNum()*(double)this->data->getNumOfChrSet());
		if(freq<this->freqthreshold)
			continue;
    	contigency[idx++]=this->Results.ControlCount[i];
    	contigency[idx++]=this->Results.CaseCount[i];
//    	std::cout<<this->Results.ControlCount[i]<<","<<this->Results.CaseCount[i]<<"\n";
    };


	  try{
		  fexact(&nrow, &validHap, contigency, &nrow, &expect, &percnt, &emin, &prt, &pre, &ws);
		  this->Results.FisherP=pre;
	  }catch(std::runtime_error &){
		  this->Results.FisherP=-1;
	  }
	  //Pearson's ChiSquare test
	  PearsonChiSquareTest(contigency,nrow,haploNum,this->Results.ChiSquare,this->Results.PearsonP);
	  delete[] contigency;
//	  std::cout<<"fisherp:"<<this->Results.FisherP;
//	  std::cout<<"\npearsonp:"<<this->Results.PearsonP;
}

std::string HaplotypeBase::reporthtml(){
	std::stringstream ss;
	ss<<"<p>Haplotypes with frequency <"<<this->freqthreshold<<" are ignored.<br>";
	ss<<"Loci chosen for haplotype analysis: ";
	for(int i=0;i<this->SnpIdx.size()-1;i++){
		ss<<this->data->vLocusName[SnpIdx[i]]<<", ";
	}
	ss<<this->data->vLocusName[SnpIdx[this->SnpIdx.size()-1]]<<"</p>\n";
	ss<<this->reporthtmltable();
	ss<<"<p><b>Global result:</b><br>Total control="<<this->data->getControlNum()<<", total case="
		<<this->data->getCaseNum()<<".<br>";
	ss<<"Global chi2 is "<<convert2string(this->Results.ChiSquare)<<", ";
	ss<<"Fisher's p is "<<convert2string(this->Results.FisherP)<<", ";
	ss<<"Pearson's p is "<<convert2string(this->Results.PearsonP)<<".</p>\n";
	return ss.str();
}

std::string HaplotypeBase::reporthtmltable(){
	boost::shared_ptr<SHEsis::CreatHtmlTable> html(new SHEsis::CreatHtmlTable());
	html->createTable("Haplotype_Analysis");
	std::vector<std::string> data;
	data.push_back("Haplotype");
	data.push_back("Case(freq)");
	data.push_back("Control(freq)");
	data.push_back("Chi2");
	data.push_back("Fisher's p");
	data.push_back("Pearson's p");
	data.push_back("OR [95% CI]");
	html->addHeadRow(data);
	for(int i=0;i<this->Results.singleHap.size();i++){
		if(-1 == this->Results.singleHap[i].pearsonp)
			continue;
		data.clear();
		std::stringstream hap;
		for(int j=0;j<this->SnpIdx.size();j++){
			hap<<this->data->getallele(this->Results.haplotypes[i][j]);
		}
		data.push_back(hap.str());
		std::string casefreq,controlfreq;
		casefreq+=convert2string((double)this->Results.CaseCount[i])+"("+
				convert2string((double)this->Results.CaseCount[i]/((double)this->data->getCaseNum()*(double)this->data->getNumOfChrSet()))+")";
		data.push_back(casefreq);
		controlfreq+=convert2string((double)this->Results.ControlCount[i])+"("+
					convert2string((double)this->Results.ControlCount[i]/((double)this->data->getControlNum()*(double)this->data->getNumOfChrSet()))+")";
		data.push_back(controlfreq);
		data.push_back(convert2string(this->Results.singleHap[i].chisquare));
		data.push_back(convert2string(this->Results.singleHap[i].pearsonp));
		data.push_back(convert2string(this->Results.singleHap[i].fisherp));
		std::string OR;
		OR+=convert2string(this->Results.singleHap[i].OR)+" ["+convert2string(this->Results.singleHap[i].orlow)
				+","+convert2string(this->Results.singleHap[i].orUp)+"]";
		data.push_back(OR);
		html->addDataRow(data);
	}
	return html->getTable();

}
}


