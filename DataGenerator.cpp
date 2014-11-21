/*
 * DataGenerator.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: ada
 */
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/shared_ptr.hpp>
#include <math.h>
#include <vector>

boost::random::mt19937 rng;
namespace po = boost::program_options;
void addOptions(int argc, char* argv[], po::options_description& desc,
                po::variables_map& vm);
void checkOptions(po::options_description& desc, po::variables_map& vm);
void printPar();
std::vector<int> GenerateQTLData(std::vector<boost::shared_ptr<short[]> > haps,
                                 int snpnum, int ploidy, int samplenum,
                                 std::string prefix);
std::vector<int> GenerateBinaryData(
    std::vector<boost::shared_ptr<short[]> > haps, int snpnum, int casenum,
    int ctrlnum, int ploidy, std::string prefix);
void printHap(std::vector<boost::shared_ptr<short[]> > haps,
              std::vector<int> count, int snpnum);
std::vector<boost::shared_ptr<short[]> > GenerateHaplotype(int alleletype,
                                                           int hapnum,
                                                           int snpnum);
struct args {
  args()
      : casenum(100),
        ctrlnum(100),
        qtlnum(0),
        ploidy(2),
        snp(10),
        output("output"),
        hap(5),
        seed(1234),
        allele(2),
        missing(0) {}
  int casenum;
  int ctrlnum;
  int qtlnum;
  int ploidy;
  int snp;
  int hap;
  int allele;
  float missing;
  std::string output;
  int seed;
} par;

int main(int argc, char* argv[]) {
  po::options_description desc("Allowed options");
  po::variables_map vm;
  addOptions(argc, argv, desc, vm);
  checkOptions(desc, vm);
  printPar();
  rng.seed(par.seed);
  std::vector<boost::shared_ptr<short[]> > haplotypes =
      GenerateHaplotype(par.allele, par.hap, par.snp);
  std::vector<int> hapcount;
  if (par.qtlnum == 0) {
    hapcount = GenerateBinaryData(haplotypes, par.snp, par.casenum, par.ctrlnum,
                                  par.ploidy, par.output);
  } else {
    hapcount = GenerateQTLData(haplotypes, par.snp, par.ploidy, par.qtlnum,
                               par.output);
  }
  printHap(haplotypes, hapcount, par.snp);
}

void addOptions(int argc, char* argv[], po::options_description& desc,
                po::variables_map& vm) {
  desc.add_options()("help", "produce help message")(
      "case", po::value<int>(), "number of cases")("ctrl", po::value<int>(),
                                                   "number of controls")(
      "qtl", po::value<int>(), "number of sample with quantitative trait")(
      "snp", po::value<int>(), "number of snps")(
      "ploidy", po::value<int>(), "number of ploidy")("hap", po::value<int>(),
                                                      "number of haplotypes")(
      "allele", po::value<int>(), "number of allele types")(
      "output", po::value<std::string>(), "output prefix")(
      "seed", po::value<int>(), "seed for random generator")(
      "missing", po::value<float>(), "missing rate");
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
}
void checkOptions(po::options_description& desc, po::variables_map& vm) {
  if (vm.count("help")) {
    std::cout << desc << "\n";
    exit(0);
  }
  if (vm.count("ploidy")) {
    par.ploidy = vm["ploidy"].as<int>();
  };
  if (vm.count("case") + vm.count("ctrl") == 1) {
    throw std::runtime_error(
        "case number and control number should both be given.");
  }
  if (vm.count("qtl") && (vm.count("case") || vm.count("ctrl"))) {
    throw std::runtime_error(
        "--qtl cannot be specified together with --case or --ctrl");
  }
  if (vm.count("case")) {
    par.casenum = vm["case"].as<int>();
  };
  if (vm.count("qtl")) {
    par.qtlnum = vm["qtl"].as<int>();
  }
  if (vm.count("ctrl")) {
    par.ctrlnum = vm["ctrl"].as<int>();
  };
  if (vm.count("snp")) {
    par.snp = vm["snp"].as<int>();
  }
  if (vm.count("seed")) {
    par.seed = vm["seed"].as<int>();
  }
  if (vm.count("hap")) {
    par.hap = vm["hap"].as<int>();
  }
  if (vm.count("allele")) {
    par.allele = vm["allele"].as<int>();
  }
  if (vm.count("output")) {
    par.output = vm["output"].as<std::string>();
  }
  if (vm.count("missing")) {
    par.missing = vm["missing"].as<float>();
  }
}

void printPar() {
  std::cout << "parameters:\n";
  if (0 == par.qtlnum) {
    std::cout << "case: " << par.casenum << ", control: " << par.ctrlnum
              << "\n";
  } else {
    std::cout << "qtl sample: " << par.qtlnum << "\n";
  }
  std::cout << "snpnum: " << par.snp << ", ploidy: " << par.ploidy
            << ", allele types: " << par.allele << "\n";
  std::cout << "seed: " << par.seed << ", hap: " << par.hap << "\n";
  std::cout << "missing rate: " << par.missing << ", output: " << par.output
            << "\n";
}

int translateRandomNumber(int number) { return (number / 10 + 1); }

std::vector<boost::shared_ptr<short[]> > GenerateHaplotype(int alleletype,
                                                           int hapnum,
                                                           int snpnum) {
  std::vector<boost::shared_ptr<short[]> > res;
  boost::random::uniform_int_distribution<> index_dist(1, alleletype * 10 - 1);
  for (int i = 0; i < hapnum; i++) {
    boost::shared_ptr<short[]> sp(new short[snpnum]);
    for (int j = 0; j < snpnum; j++) {
      sp[j] = (short)translateRandomNumber(index_dist(rng));
    }  //
    res.push_back(sp);
  }
  return res;
}

std::vector<int> GenerateBinaryData(
    std::vector<boost::shared_ptr<short[]> > haps, int snpnum, int casenum,
    int ctrlnum, int ploidy, std::string prefix) {
  std::vector<int> hapcount;
  hapcount.resize(haps.size(), 0);
  std::ofstream casefile;
  std::ofstream ctrlfile;
  casefile.open((prefix + "_case.txt").c_str());
  if (!casefile)
    throw std::runtime_error("Unable to open output file: " + prefix +
                             "_case.txt");
  ctrlfile.open((prefix + "_ctrl.txt").c_str());
  if (!ctrlfile)
    throw std::runtime_error("Unable to open output file: " + prefix +
                             "_ctrl.txt");
  std::vector<int> selection;
  boost::random::uniform_int_distribution<> HapSelection(0, haps.size() - 1);
  boost::random::uniform_int_distribution<> Missing(0, 1000);
  for (int i = 0; i < casenum; i++) {
    casefile << "case" << i << " ";
    selection.resize(ploidy, -1);
    for (int p = 0; p < ploidy; p++) {
      int idx = HapSelection(rng);
      selection[p] = idx;
      hapcount[idx]++;
    }
    for (int s = 0; s < snpnum; s++) {
      std::vector<int> tmp;
      for (int select = 0; select < selection.size(); select++) {
        BOOST_ASSERT(selection[select] != -1);
        if (Missing(rng) < par.missing * 1000)
          tmp.push_back(0);
        else
          tmp.push_back(haps[selection[select]][s]);
      }
      std::random_shuffle(tmp.begin(), tmp.end());
      for (int select = 0; select < selection.size(); select++) {
        casefile << tmp[select] << " ";
      }
    }
    casefile << "\n";
  }
  casefile.close();

  for (int i = 0; i < ctrlnum; i++) {
    ctrlfile << "ctrl" << i << " ";
    selection.resize(ploidy, -1);
    for (int p = 0; p < ploidy; p++) {
      int idx = HapSelection(rng);
      selection[p] = idx;
      hapcount[idx]++;
    }
    for (int s = 0; s < snpnum; s++) {
      std::vector<int> tmp;
      for (int select = 0; select < selection.size(); select++) {
        BOOST_ASSERT(selection[select] != -1);
        if (Missing(rng) < par.missing * 1000)
          tmp.push_back(0);
        else
          tmp.push_back(haps[selection[select]][s]);
      }
      std::random_shuffle(tmp.begin(), tmp.end());
      for (int select = 0; select < selection.size(); select++) {
        ctrlfile << tmp[select] << " ";
      }
    }
    ctrlfile << "\n";
  }
  ctrlfile.close();
  return hapcount;
}

std::vector<int> GenerateQTLData(std::vector<boost::shared_ptr<short[]> > haps,
                                 int snpnum, int ploidy, int samplenum,
                                 std::string prefix) {
  std::vector<int> hapcount;
  hapcount.resize(haps.size(), 0);
  std::ofstream qtlfile;
  qtlfile.open((prefix + "_qtl.txt").c_str());
  if (!qtlfile)
    throw std::runtime_error("Unable to open output file: " + prefix +
                             "_case.txt");
  std::vector<int> selection;
  boost::random::uniform_int_distribution<> HapSelection(0, haps.size() - 1);
  boost::random::uniform_int_distribution<> Missing(0, 1000);
  boost::normal_distribution<> normal(1, 10);
  for (int i = 0; i < samplenum; i++) {
    qtlfile << "qtl" << i << " " << abs(normal(rng)) << " ";
    selection.resize(ploidy, -1);
    for (int p = 0; p < ploidy; p++) {
      int idx = HapSelection(rng);
      selection[p] = idx;
      hapcount[idx]++;
    }
    for (int s = 0; s < snpnum; s++) {
      std::vector<int> tmp;
      for (int select = 0; select < selection.size(); select++) {
        BOOST_ASSERT(selection[select] != -1);
        if (Missing(rng) < par.missing * 1000)
          tmp.push_back(0);
        else
          tmp.push_back(haps[selection[select]][s]);
      }
      std::random_shuffle(tmp.begin(), tmp.end());
      for (int select = 0; select < selection.size(); select++) {
        qtlfile << tmp[select] << " ";
      }
    }
    qtlfile << "\n";
  }
  qtlfile.close();
  return hapcount;
}

void printHap(std::vector<boost::shared_ptr<short[]> > haps,
              std::vector<int> count, int snpnum) {
  for (int hap = 0; hap < haps.size(); hap++) {
    for (int snp = 0; snp < snpnum; snp++) {
      std::cout << haps[hap][snp];
    }
    std::cout << " :" << count[hap] << "\n";
  }
}
