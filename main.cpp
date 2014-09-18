/*
 * main.cpp
 *
 *  Created on: Aug 3, 2014
 *      Author: ionadmin
 */
#include "SHEsisData.h"
#include "Haplotype.h"
#include "HaplotypeDiploid.h"
#include "AssociationTest.h"
#include "LDTest.h"
#include "HWETest.h"
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <fstream>

namespace po = boost::program_options;

typedef enum{
	LD_IN_CASE,
	LD_IN_CTRL,
	LD_IN_BOTH
}LD_TYPE;

struct arguments{
	arguments():haploAnalysis(false),assocAnalysis(false),hweAnalysis(false),ldAnalysis(false),permutation(-1){};
	std::vector<std::string> inputfiles;
	std::vector<std::string> inputcases;
	std::vector<std::string> inputctrls;
	int ploidy;
	int permutation;
	bool containsPhenotype;
	bool haploAnalysis;
	bool assocAnalysis;
	bool hweAnalysis;
	bool ldAnalysis;
	std::vector<short> mask;
	LD_TYPE ldtype;
} SHEsisArgs;


void addOptions(int argc, char *argv[],po::options_description& desc,po::variables_map& vm );
void checkOptions(po::options_description& desc,po::variables_map& vm);
int ReadInput(int ploidy, bool containsPhenotype,std::string filepath, std::vector<std::vector<std::string> >& filecontent);
boost::shared_ptr<SHEsis::SHEsisData> parseDataWithPhenotype(int snpnum,int ploidy, std::vector<std::vector<std::string> >& content);
boost::shared_ptr<SHEsis::SHEsisData> parseDataNoPhenotype(int snpnum,int ploidy,
		std::vector<std::vector<std::string> >& casecontent,
		std::vector<std::vector<std::string> >& ctrlcontent);
boost::shared_ptr<SHEsis::SHEsisData> parseInput();


int main(int argc, char *argv[])
{
	po::options_description desc("Allowed options");
	po::variables_map vm;
	boost::shared_ptr<SHEsis::SHEsisData> data;
	addOptions(argc,argv,desc,vm);
	try{
		checkOptions(desc,vm);
		data=parseInput();
	}catch(std::runtime_error& e){
		std::cout<<"***ERROR:"<<e.what()<<"\n";
		std::cout<<desc<<"\n";
		exit(-1);
	};
	boost::shared_ptr<SHEsis::AssociationTest> AssocHandle;//(new SHEsis::AssociationTest(data));
	boost::shared_ptr<SHEsis::HWETest> HWEHandle;//(new SHEsis::HWETest(data));
	boost::shared_ptr<SHEsis::HaplotypeBase> HapHandle;//(new SHEsis::Haplotype(data));
	//boost::shared_ptr<SHEsis::HaplotypeDiploid> DiploidHapHandle;//(new SHEsis::Haplotype(data));
	boost::shared_ptr<SHEsis::LDTest> LDHandle;

	if(SHEsisArgs.assocAnalysis){
		AssocHandle.reset(new SHEsis::AssociationTest(data));
		if(SHEsisArgs.permutation!=-1){
			AssocHandle->setPermutationTimes(SHEsisArgs.permutation);
			AssocHandle->permutation();
		}else{
			AssocHandle->association();
		}
	};

	if(SHEsisArgs.hweAnalysis){
		HWEHandle.reset(new SHEsis::HWETest(data));
		HWEHandle->AllSnpHWETest();
	}

	if(SHEsisArgs.haploAnalysis){
		if(data->getNumOfChrSet()<=2){
			if(SHEsisArgs.mask.size()!=data->getSnpNum()){
				HapHandle.reset(new SHEsis::HaplotypeDiploid(data));
			}else{
				int snpnum=0;
				for(int i=0;i<SHEsisArgs.mask.size();i++){
					snpnum+=SHEsisArgs.mask[i];
				}
				HapHandle.reset(new SHEsis::HaplotypeDiploid(data,snpnum,SHEsisArgs.mask));
			}
		}else{
			if(SHEsisArgs.mask.size()!=data->getSnpNum()){
				HapHandle.reset(new SHEsis::Haplotype(data));
			}else{
				int snpnum=0;
				for(int i=0;i<SHEsisArgs.mask.size();i++){
					snpnum+=SHEsisArgs.mask[i];
				}
				HapHandle.reset(new SHEsis::Haplotype(data,snpnum,SHEsisArgs.mask));
			}
		}
		HapHandle->startHaplotypeAnalysis();
	};

	if(SHEsisArgs.ldAnalysis){
		LDHandle.reset(new SHEsis::LDTest(data));
		LDHandle->AllLociLDtest();
	};

	return 0;
}

boost::shared_ptr<SHEsis::SHEsisData> parseInput(){
	boost::shared_ptr<SHEsis::SHEsisData> data;
	std::vector<std::vector<std::string> > filecontent;
	std::vector<std::vector<std::string> > filecontentcase;
	std::vector<std::vector<std::string> > filecontentctrl;
	int snpnum,snpcase,snpctrl;
	if(SHEsisArgs.inputfiles.size()>0){
		for(int i=0;i<SHEsisArgs.inputfiles.size();i++){
			snpnum=ReadInput(SHEsisArgs.ploidy,SHEsisArgs.containsPhenotype,SHEsisArgs.inputfiles[i],filecontent);
		}
		data=parseDataWithPhenotype(snpnum,SHEsisArgs.ploidy,filecontent);
	}else if(SHEsisArgs.inputctrls.size()>0 && SHEsisArgs.inputcases.size()>0){
		for(int i=0;i<SHEsisArgs.inputctrls.size();i++)
			snpctrl=ReadInput(SHEsisArgs.ploidy,SHEsisArgs.containsPhenotype,SHEsisArgs.inputctrls[i],filecontentctrl);
		for(int i=0;i<SHEsisArgs.inputcases.size();i++)
			snpcase=ReadInput(SHEsisArgs.ploidy,SHEsisArgs.containsPhenotype,SHEsisArgs.inputcases[i],filecontentcase);
		if(snpctrl != snpcase)
			throw std::runtime_error("SNP number in cases and controls disagrees.");
		data=parseDataNoPhenotype(snpcase,SHEsisArgs.ploidy,filecontentcase,filecontentctrl);
	}else{
		throw std::runtime_error("Error:Please check the input files.");
	};
	return data;
}

boost::shared_ptr<SHEsis::SHEsisData> parseDataWithPhenotype(int snpnum,int ploidy, std::vector<std::vector<std::string> >& content){
	int samplenum=content.size();
	boost::shared_ptr<SHEsis::SHEsisData> pdata(new SHEsis::SHEsisData(samplenum,snpnum,ploidy));
	for(int sample=0;sample<samplenum;sample++){
		pdata->vLabel[sample]=(SHEsis::SampleStatus)std::atoi(content[sample][1].c_str());
		for(int snp=0;snp<snpnum;snp++){
			for(int p=0;p<ploidy;p++){
				pdata->mGenotype[sample][snp][p]=pdata->GetAlleleCode(content[sample][snp*ploidy+2+p]);
			}
		}
	};
	return pdata;
}

boost::shared_ptr<SHEsis::SHEsisData> parseDataNoPhenotype(int snpnum,int ploidy,
		std::vector<std::vector<std::string> >& casecontent,
		std::vector<std::vector<std::string> >& ctrlcontent){
	int casenum=casecontent.size();
	int ctrlnum=ctrlcontent.size();
	boost::shared_ptr<SHEsis::SHEsisData> pdata(new SHEsis::SHEsisData(casenum+ctrlnum,snpnum,ploidy));
	for(int sample=0;sample<casecontent.size();sample++){
		pdata->vLabel[sample]=SHEsis::CASE;
		for(int snp=0;snp<snpnum;snp++){
			for(int p=0;p<ploidy;p++){
				pdata->mGenotype[sample][snp][p]=pdata->GetAlleleCode(casecontent[sample][snp*ploidy+1+p]);
			}
		}
	};
	for(int sample=0;sample<casecontent.size();sample++){
		pdata->vLabel[sample+casenum]=SHEsis::CONTROL;
		for(int snp=0;snp<snpnum;snp++){
			for(int p=0;p<ploidy;p++){
				pdata->mGenotype[sample+casenum][snp][p]=pdata->GetAlleleCode(casecontent[sample][snp*ploidy+1+p]);
			}
		}
	};
	return pdata;
}

void addOptions(int argc, char *argv[],po::options_description& desc,po::variables_map& vm ){
	desc.add_options()
	    ("help", "produce help message")
	    ("input",po::value<std::vector<std::string> >(),"path for the input file containing both cases and controls")
	    ("input-case",po::value<std::vector<std::string> >(),"path for the input file containing cases")
	    ("input-ctrl",po::value<std::vector<std::string> >(),"path for the input file containing controls")
	    ("ploidy",po::value<int>(),"number of ploidy")
	    ("assoc","perform case/control association test")
	    ("permutation",po::value<int>(),"times for permutation")
	    ("haplo","perform haplotype analysis")
	    //("error-rate",po::value<double>()->default_value(0.00001),"convergence threshold for haplotype analysis")
	    ("mask",po::value<std::string>(),"mask of snps for haplotype analysis, comma delimited. eg. mask=1,0,1 to use 1st and 3rd SNPs when there are 3 SNPs in all.")
	    ("ld-in-case","perform Linkage disequilibrium test in cases")
	    ("ld-in-ctrl","perform Linkage disequilibrium test in controls")
	    ("ld","perform Linkage disequilibrium test in both cases and controls")
	    ("hwe","perform Hardy-Weinberg disequilibrium test")
	    ;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
}

void checkOptions(po::options_description& desc,po::variables_map& vm){
	if(vm.count("help")){
		std::cout<<desc<<"\n";
		exit(0);
	}

	if((0==vm.count("input"))&&(0==vm.count("input-case"))&&(0==vm.count("input-ctrl")))
		throw std::runtime_error("no input file specified.");
	if(vm.count("input")&&(vm.count("input-case")||vm.count("input-ctrl")))
		throw std::runtime_error("--input and --input-case/--input-ctrl cannot be specified at the same time.");
	if((0==vm.count("input-case")&&(0!=vm.count("input-ctrl"))) || (0==vm.count("input-ctrl")&&(0!=vm.count("input-case"))))
		throw std::runtime_error("-input-case and --input-ctrl should both be specified.");
	if(vm.count("input")>0){
		SHEsisArgs.inputfiles=vm["input"].as< std::vector<std::string> >();
		SHEsisArgs.containsPhenotype=true;
	}else if(vm.count("input-case")>0 && vm.count("input-ctrl")>0){
		SHEsisArgs.inputcases=vm["input-case"].as< std::vector<std::string> >();
		SHEsisArgs.inputctrls=vm["input-ctrl"].as< std::vector<std::string> >();
		SHEsisArgs.containsPhenotype=false;
	}else{
		throw std::runtime_error("error in parsing input files.");
	}

	if(0 == vm.count("ploidy"))
		throw std::runtime_error("no ploidy information given.");
	else if(vm["ploidy"].as<int>() <=0)
		throw std::runtime_error("number of ploidy should be higher than 0.");
	else
		SHEsisArgs.ploidy=vm["ploidy"].as<int>();

	int ldcase=vm.count("ld-in-case")>0?1:0;
	int ldctrl=vm.count("ld-in-ctrl")>0?1:0;
	int ld=vm.count("ld")>0?1:0;
	if(ldcase+ldctrl+ld>1){
		throw std::runtime_error("--ld-in-case/--ld-in-ctrl/--ld cannot be specified at the same time");
	}else if(ldcase+ldctrl+ld == 1)
	{
		SHEsisArgs.ldAnalysis=true;
		if(ldcase)
			SHEsisArgs.ldtype=LD_IN_CASE;
		if(ldctrl)
			SHEsisArgs.ldtype=LD_IN_CTRL;
		if(ld)
			SHEsisArgs.ldtype=LD_IN_BOTH;
	};

	if(vm.count("permutation")>0 && vm.count("assoc") == 0)
		throw std::runtime_error("--permutaion should be used along with --assoc");
	if(vm.count("assoc")){
		SHEsisArgs.assocAnalysis=true;
		if(vm.count("permutation")){
			SHEsisArgs.permutation=vm["permutation"].as<int>();
			if(SHEsisArgs.permutation<=0)
				throw std::runtime_error("permutation times should be higher than 0.");
		}
	};

	if(vm.count("hwe")){
		SHEsisArgs.hweAnalysis=true;
	};
	if(vm.count("haplo")){
		SHEsisArgs.haploAnalysis=true;
	}

	if(vm.count("mask")){
		std::string maskstr=vm["mask"].as<std::string>();
		std::vector<std::string> maskvec;
		boost::split(maskvec,maskstr,boost::is_any_of("\t ,"));
		for(int i=0;i<maskvec.size();i++){
			if(std::strcmp("0",maskvec[i].c_str()) ==0 ){
				SHEsisArgs.mask.push_back(0);
			}else if(std::strcmp("1",maskvec[i].c_str()) ==0){
				SHEsisArgs.mask.push_back(1);
			}else{
				throw std::runtime_error("mask should be composed of 0 and 1. But "+maskvec[i]+" found.");
			}
		}
	}

	if(!SHEsisArgs.hweAnalysis && !SHEsisArgs.haploAnalysis && !SHEsisArgs.assocAnalysis && !SHEsisArgs.ldAnalysis)
		throw std::runtime_error("at least one type of analysis should be specified.");
}


int ReadInput(int ploidy, bool containsPhenotype,std::string filepath, std::vector<std::vector<std::string> >& filecontent){
	std::string line;
	std::ifstream file(filepath);
	if(!file.is_open())
		throw std::runtime_error("Cannot open input file.");
	int lineidx=1;
	int snpnum=0;
	int expectedfileds=0;
	while(getline(file,line)){
		if(line.empty()){
			lineidx++;
			continue;
		}
		std::vector<std::string> strs;
		boost::split(strs,line,boost::is_any_of("\t ,"));

		if(containsPhenotype){
			if((SHEsis::SampleStatus)(std::atoi(strs[1].c_str())!=SHEsis::CASE &&
				(SHEsis::SampleStatus)(std::atoi(strs[1].c_str())) !=SHEsis::CONTROL)){
				std::stringstream ss;
				ss<<"Error in line "<<lineidx<<" ,phenotype should be either "<<SHEsis::CASE<<" or "<<SHEsis::CONTROL<<", but "<<strs[1]<<" found";
				throw std::runtime_error(ss.str());
			};
		};

		if(lineidx == 1){
			if((strs.size()-1-(int)containsPhenotype)%ploidy != 0)
				throw std::runtime_error("Error in line 1, file:"+filepath+". Either snp num or ploidy num is wrong");
			snpnum=(strs.size()-1-(int)containsPhenotype)/ploidy;
			expectedfileds=(1+(int)containsPhenotype+snpnum*ploidy);
			if(filecontent.size()!=0 && filecontent[0].size()!=expectedfileds){
				throw std::runtime_error("number of snps disagree between files");
			}
		}

		if(strs.size()!=expectedfileds){
			std::stringstream ss;
			ss<<"Error in line "<<lineidx<<", file"<<filepath<<", expecting "<<expectedfileds<<" fileds, but"<<strs.size()<<" fileds are found";
			throw std::runtime_error(ss.str());
		};

		filecontent.push_back(strs);
	}

	return snpnum;
}




