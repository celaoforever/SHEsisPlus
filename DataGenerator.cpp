/*
 * DataGenerator.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: ada
 */
#include <boost/program_options.hpp>
#include <fstream>
#include <string>
namespace po = boost::program_options;
void addOptions(int argc, char* argv[], po::options_description& desc,
                po::variables_map& vm);
void checkOptions(po::options_description& desc, po::variables_map& vm);
struct args{
	args():casenum(100),ctrlnum(100),qtlnum(0),ploidy(2),snp(10),output("output"),seed(1234){}
	int casenum;
	int ctrlnum;
	int qtlnum;
	int ploidy;
	int snp;
	std::string output;
	int seed;
}par;

int main(int argc, char* argv[]) {
  po::options_description desc("Allowed options");
  po::variables_map vm;
  boost::shared_ptr<SHEsis::SHEsisData> data;
  addOptions(argc, argv, desc, vm);
}

void addOptions(int argc, char* argv[], po::options_description& desc,
                po::variables_map& vm){
	  desc.add_options()("help", "produce help message")
						("case", po::value<int>(), "number of cases")
						("ctrl", po::value<int>(), "number of controls")
						("qtl",  po::value<int>(), "number of sample with quantitative trait")
						("snp",  po::value<int>(), "number of snps")
						("ploidy",po::value<int>(),"number of ploidy")
						("output",po::value<std::string>(),"output prefix")
						("seed", po::value<int>(),"seed for random generator");
	  po::store(po::parse_command_line(argc, argv, desc), vm);
	  po::notify(vm);
}
void checkOptions(po::options_description& desc, po::variables_map& vm){
	  if (vm.count("help")) {
	    std::cout << desc << "\n";
	    exit(0);
	  }
	  if(vm.count("ploidy")){
		  par.ploidy=vm["ploidy"].as<int>();
	  };
	  if(vm.count("case")+vm.count("ctrl") == 1){
		  throw std::runtime_error("case number and control number should both be given.");
	  }
	  if(vm.count("qtl") && (vm.count("case")||vm.count("ctrl"))){
		  throw std::runtime_error("--qtl cannot be specified together with --case or --ctrl");
	  }
	  if(vm.count("case")){
		  par.casenum=vm["case"].as<int>();
	  };
	  if(vm.count("ctrl")){
		  par.ctrlnum=vm["ctrl"].as<int>();
	  };
	  if(vm.count("snp")){
		  par.snp=vm["snp"].as<int>();
	  }
	  if(vm.count("seed")){
		  par.seed=vm["seed"].as<int>();
	  }
	  if(vm.count("output")){
		  par.output=vm["output"].as<std::string>();
	  }
}

void printPar(){

}
