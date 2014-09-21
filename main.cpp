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
#ifdef _WIN32
	#include <windows.h>
#endif
#define VERSION 2.0

namespace po = boost::program_options;
std::stringstream report;

struct arguments{
	arguments():haploAnalysis(false),assocAnalysis(false),html(true),hweAnalysis(false),ldAnalysis(false),permutation(-1),lft(0.03){};
	std::vector<std::string> inputfiles;
	std::vector<std::string> inputcases;
	std::vector<std::string> inputctrls;
	std::vector<std::string> snpnames;
	std::string output;
	int ploidy;
	int permutation;
	bool containsPhenotype;
	bool haploAnalysis;
	bool assocAnalysis;
	bool hweAnalysis;
	bool ldAnalysis;
	double lft;
	std::vector<short> mask;
	SHEsis::LD_TYPE ldtype;
	bool html;
} SHEsisArgs;


void addOptions(int argc, char *argv[],po::options_description& desc,po::variables_map& vm );
void checkOptions(po::options_description& desc,po::variables_map& vm);
int ReadInput(int ploidy, bool containsPhenotype,std::string filepath, std::vector<std::vector<std::string> >& filecontent);
boost::shared_ptr<SHEsis::SHEsisData> parseDataWithPhenotype(int snpnum,int ploidy, std::vector<std::vector<std::string> >& content);
boost::shared_ptr<SHEsis::SHEsisData> parseDataNoPhenotype(int snpnum,int ploidy,
		std::vector<std::vector<std::string> >& casecontent,
		std::vector<std::vector<std::string> >& ctrlcontent);
boost::shared_ptr<SHEsis::SHEsisData> parseInput();
void getSnpNamefile(std::string& path,std::vector<std::string>& snp);
void getSnpNameline(std::string& names,std::vector<std::string>& snp);
void writeHtmlHeader();
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
		std::cout<<"***ERROR: "<<e.what()<<"\n";
		std::cout<<desc<<"\n";
		exit(-1);
	};
	if((SHEsisArgs.snpnames.size()!=0 )&& (data->getSnpNum() !=SHEsisArgs.snpnames.size())){
		std::cout<<"***WARNING: given number of snp names is not the same at that in data. Ignoring snp names...\n";
		SHEsisArgs.snpnames.clear();
	};
	for(int i=0;i<SHEsisArgs.snpnames.size();i++){
		data->setLocusName(i,SHEsisArgs.snpnames[i]);
	}
	if(SHEsisArgs.html){
		writeHtmlHeader();
	}

	std::ofstream ofile;
	std::string filename=SHEsisArgs.output+(SHEsisArgs.html?".html":".txt");
	ofile.open(filename.c_str());
	if(!ofile){
		std::cout<<"***ERROR: cannot open output file "<<filename<<"\n";
		exit(-1);
	}

	boost::shared_ptr<SHEsis::AssociationTest> AssocHandle;
	boost::shared_ptr<SHEsis::HWETest> HWEHandle;
	boost::shared_ptr<SHEsis::HaplotypeBase> HapHandle;
	boost::shared_ptr<SHEsis::LDTest> LDHandle;

	report<<"<h1>SHEsis </h1>\n";

	if(SHEsisArgs.assocAnalysis){
		std::cout<<"Starting association test...\n";
		AssocHandle.reset(new SHEsis::AssociationTest(data));
		if(SHEsisArgs.permutation!=-1){
			AssocHandle->setPermutationTimes(SHEsisArgs.permutation);
			AssocHandle->permutation();
		}else{
			AssocHandle->association();
		}
		report<<"\n<h2> Single Locus Association Test: </h2>\n";
		report<<AssocHandle->reporthtml();
		std::cout<<"done\n";
	};

	if(SHEsisArgs.hweAnalysis){
		std::cout<<"Starting Hardy-Weinberg equilibrium test...\n";
		HWEHandle.reset(new SHEsis::HWETest(data));
		HWEHandle->AllSnpHWETest();
		report<<"\n<h2> Hardy-Weinberg Equilibrium Test: </h2>\n";
		report<<HWEHandle->reporthtml(0.1);
		std::cout<<"done\n";
	}

	if(SHEsisArgs.haploAnalysis){
		std::cout<<"Starting haplotype analysis...\n";
		if(data->getNumOfChrSet()<=2){
			if(SHEsisArgs.mask.size()!=data->getSnpNum()){
				if(SHEsisArgs.mask.size()!=0)
					std::cout<<"***WARNING: length of given mask is wrong. Ignoring mask...\n";
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
				if(SHEsisArgs.mask.size()!=0)
					std::cout<<"***WARNING: length of given mask is wrong. Ignoring mask...\n";
				HapHandle.reset(new SHEsis::Haplotype(data));
			}else{
				int snpnum=0;
				for(int i=0;i<SHEsisArgs.mask.size();i++){
					snpnum+=SHEsisArgs.mask[i];
				}
				HapHandle.reset(new SHEsis::Haplotype(data,snpnum,SHEsisArgs.mask));
			}
		}
		if(SHEsisArgs.lft>0 && SHEsisArgs.lft <1)
			HapHandle->setFreqThreshold(SHEsisArgs.lft);
		else
			std::cout<<"***WARNING: lowest frequency threshold for haplotype analysis is invalid..defaulting to 0.03\n";
		HapHandle->setSilent(false);
		HapHandle->startHaplotypeAnalysis();
		HapHandle->AssociationTest();
		report<<"\n<h2> Haplotype Analysis: </h2>\n";
		report<<HapHandle->reporthtml();
		std::cout<<"done\n";
	};

	if(SHEsisArgs.ldAnalysis){
		std::cout<<"Starting linkage disequilibrium analysis...\n";
		LDHandle.reset(new SHEsis::LDTest(data,SHEsisArgs.output+".bmp"));
		LDHandle->AllLociLDtest();
		LDHandle->DrawLDMap();
		report<<"\n<h2> Linkage Disequilibrium Analysis: </h2>\n";
		report<<LDHandle->reporthtml();
		std::cout<<"done\n";
	};
	if(SHEsisArgs.html){
		report<<"</body>\n</html>\n";
	};

	ofile<<report.str();
	ofile.close();
	std::string cmd="firefox "+SHEsisArgs.output+".html";
#ifdef _WIN32
	ShellExecute(NULL,"open",SHEsisArgs.output+".html",NULL,NULL,SW_SHOWNORMAL);
#else
	system(cmd.c_str());
#endif

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
		throw std::runtime_error("Please check the input files.");
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
	for(int sample=0;sample<ctrlcontent.size();sample++){
		pdata->vLabel[sample+casenum]=SHEsis::CONTROL;
		for(int snp=0;snp<snpnum;snp++){
			for(int p=0;p<ploidy;p++){
				pdata->mGenotype[sample+casenum][snp][p]=pdata->GetAlleleCode(ctrlcontent[sample][snp*ploidy+1+p]);
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
	    ("snpname-file",po::value<std::string>(),"path for file that contains names of snps")
	    ("snpname-line",po::value<std::string>(),"snp names are as arguments")
	    ("output",po::value<std::string>(),"prefix of output files")
	    ("ploidy",po::value<int>(),"number of ploidy")
	    ("assoc","perform case/control association test")
	    ("permutation",po::value<int>(),"times for permutation")
	    ("haplo","perform haplotype analysis")
	    //("error-rate",po::value<double>()->default_value(0.00001),"convergence threshold for haplotype analysis")
	    ("mask",po::value<std::string>(),"mask of snps for haplotype analysis, comma delimited. eg. mask=1,0,1 to use 1st and 3rd SNPs when there are 3 SNPs in all.")
	    ("lft",po::value<double>(),"lowest frequency threshold for haplotype analysis")
	    ("ld-in-case","perform Linkage disequilibrium test in cases")
	    ("ld-in-ctrl","perform Linkage disequilibrium test in controls")
	    ("ld","perform Linkage disequilibrium test in both cases and controls")
	    ("hwe","perform Hardy-Weinberg disequilibrium test")
	    ;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
}

void getSnpNamefile(std::string path,std::vector<std::string>& snp){
	std::string line;
	std::vector<std::string> strs;
	std::ifstream file(path.c_str());
	if(!file.is_open())
		throw std::runtime_error("Cannot open snp name file.");
	while(getline(file,line)){
		if(line.empty())
			continue;
		boost::trim_if(line,boost::is_any_of("\t "));
		boost::split(strs,line,boost::is_any_of("\t "),boost::token_compress_on);
		if(strs.size()>1){
			std::cout<<"***WARNING: "<<path<<" contains more than 1 fileds in a single line.";
			std::cout<<"Ignoring this file...\n";
			snp.clear();
			return;
		}
		snp.push_back(line);
	}
}

void getSnpNameline(std::string names,std::vector<std::string>& snp){
	boost::trim_if(names,boost::is_any_of("\t ,"));
	boost::split(snp,names,boost::is_any_of("\t ,"),boost::token_compress_on);
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
		throw std::runtime_error("--input-case and --input-ctrl should both be specified.");
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
	if(0!=vm.count("snpname-file") && 0!=vm.count("snpname-line"))
		throw std::runtime_error("--snpname-line and --snpname-file cannot be specified at the same time.");
	if(vm.count("snpname-file")){
		//SHEsisArgs.snppath=vm["snpname-file"].as<std::string>();
		getSnpNamefile(vm["snpname-file"].as<std::string>(),SHEsisArgs.snpnames);
	}else if(vm.count("snpname-line")){
		//SHEsisArgs.snpline=vm["snpname-line"].as<std::string>();
		getSnpNameline(vm["snpname-line"].as<std::string>(),SHEsisArgs.snpnames);
	}

	if(0==vm.count("output"))
		SHEsisArgs.output="SHEsis";
	else
		SHEsisArgs.output=vm["output"].as<std::string>();

	if(0!=vm.count("lft"))
		SHEsisArgs.lft=vm["lft"].as<double>();

	if(0 == vm.count("ploidy"))
		throw std::runtime_error("no ploidy information given.");
	else if(vm["ploidy"].as<int>() <=0)
		throw std::runtime_error("number of ploidy should be higher than 0.");
	else
		SHEsisArgs.ploidy=vm["ploidy"].as<int>();

	int ldcase=vm.count("ld-in-case");
	int ldctrl=vm.count("ld-in-ctrl");
	int ld=vm.count("ld");
	if(ldcase+ldctrl+ld>1){
		throw std::runtime_error("--ld-in-case/--ld-in-ctrl/--ld cannot be specified at the same time");
	}else if(ldcase+ldctrl+ld == 1)
	{
		SHEsisArgs.ldAnalysis=true;
		if(ldcase)
			SHEsisArgs.ldtype=SHEsis::LD_IN_CASE;
		if(ldctrl)
			SHEsisArgs.ldtype=SHEsis::LD_IN_CTRL;
		if(ld)
			SHEsisArgs.ldtype=SHEsis::LD_IN_BOTH;
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
	std::ifstream file(filepath.c_str());
	if(!file.is_open())
		throw std::runtime_error("Cannot open input file: "+filepath);
	int lineidx=1;
	int snpnum=0;
	int expectedfileds=0;
	while(getline(file,line)){
		if(line.empty()){
			lineidx++;
			continue;
		}
		std::vector<std::string> strs;
		boost::trim_if(line,boost::is_any_of("\t ,"));
		boost::split(strs,line,boost::is_any_of("\t ,"),boost::token_compress_on);

		if(containsPhenotype){
			if((SHEsis::SampleStatus)(std::atoi(strs[1].c_str())!=SHEsis::CASE &&
				(SHEsis::SampleStatus)(std::atoi(strs[1].c_str())) !=SHEsis::CONTROL)){
				std::stringstream ss;
				ss<<"Error in line "<<lineidx<<", file:"<<filepath<<", phenotype should be either "<<SHEsis::CASE<<" for cases or "<<SHEsis::CONTROL<<" for controls, but "<<strs[1]<<" found";
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
			ss<<"Error in line "<<lineidx<<", file:"<<filepath<<", expecting "<<expectedfileds<<" fileds, but "<<strs.size()<<" fileds are found";
			throw std::runtime_error(ss.str());
		};
		lineidx++;
		filecontent.push_back(strs);
	}

	return snpnum;
}

void writeHtmlHeader(){
	report<<"<!DOCTYPE html> \n\
	<html> \n\
	<head> \n\
	<style> \n\
	table { \n\
	    width:auto; \n\
	} \n\
	table, th, td { \n\
	    border: 1px solid black; \n\
	    border-collapse: collapse; \n\
	} \n\
	th, td { \n\
	    padding: 5px; \n\
	    text-align: left; \n\
	} \n\
	table tr:nth-child(even) { \n\
	    background-color: #eee;\n\
	}\n\
	table tr:nth-child(odd) {\n\
	   background-color:#fff;\n\
	}\n\
	table th    {\n\
	    background-color: #AAA;\n\
	    color: white;\n\
	}\n\
	</style>\n\
	</head>\n\
\n\
	<body>\n";
}





