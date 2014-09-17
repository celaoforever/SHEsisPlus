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

void parseOptions(int argc, char *argv[],po::options_description& desc,po::variables_map& vm );

int main(int argc, char *argv[])
{
	po::options_description desc("Allowed options");
	po::variables_map vm;
	parseOptions(argc,argv,desc,vm);



	return 0;
}

void parseOptions(int argc, char *argv[],po::options_description& desc,po::variables_map& vm ){
	desc.add_options()
	    ("help", "produce help message")
	    ("input",po::value<std::string>(),"path for the input file containing both cases and controls")
	    ("input-case",po::value<std::string>(),"path for the input file containing cases")
	    ("input-ctrl",po::value<std::string>(),"path for the input file containing controls")
	    ("ploidy",po::value<int>(),"number of ploidy")
	    ("assoc","perform case/control association test")
	    ("haplo","perform haplotype analysis")
	    ("mask",po::value<std::string>(),"mask of snps for haplotype analysis, eg. mask=101 to use 1st and 3rd SNPs.")
	    ("ld",po::value<std::string>(),"perform Linkage disequilibrium test in case/control/both, allowed value:{CASE,CONTROL,BOTH},default:BOTH")
	    ("hwe",po::value<std::string>(),"perform Hardy-weinberg disequilibrium test in case/control/both, allowed value:{CASE,CONTROL,BOTH},default:BOTH")
	    ;
}

void ReadInput(int ploidy, bool containsPhenotype,std::string filepath, boost::shared_ptr<SHEsis::SHEsisData>& data){
	std::string line;
	std::ifstream file(filepath);
	std::vector<std::string> filecontent;
	if(!file.is_open())
		throw std::runtime_error("Cannot open input file.");
	int lineidx=1;
	//int sampleidx=0;
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
			if(strs[1]!=SHEsis::CASE && strs[1] !=SHEsis::CONTROL){
				std::stringstream ss;
				ss<<"Error in line "<<lineidx<<" ,phenotype should be either "<<SHEsis::CASE<<" or "<<SHEsis::CONTROL<<", but "<<strs[1]<<" found";
				throw std::runtime_error(ss.str());
			};
		};

		if(lineidx == 1){
			if((strs.size()-1-(int)containsPhenotype)%ploidy != 0)
				throw std::runtime_error("Error in line 1. Either snp num or ploidy num is wrong");
			snpnum=(strs.size()-1-(int)containsPhenotype)/ploidy;
			expectedfileds=(1+(int)containsPhenotype+snpnum*ploidy);
		}

		if(strs.size()!=expectedfileds){
			std::stringstream ss;
			ss<<"Error in line "<<lineidx<<", expecting "<<expectedfileds<<" fileds, but"<<strs.size()<<" fileds are found";
			throw std::runtime_error(ss.str());
		};

		filecontent.push_back(line);
	}

//	data.reset(new SHEsis::SHEsisData(filecontent.size(),snpnum,ploidy));
//	for(int i=0;i<filecontent.size())


}




