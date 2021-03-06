/*
 * main.cpp
 *
 *  Created on: Aug 3, 2014
 *      Author: ionadmin
 */
#include "SHEsisData.h"
#include "HaplotypeEM.h"
#include "Haplotype.h"
#include "HaplotypeLD.h"
//#include "HaplotypeDiploid.h"
#include "AssociationTest.h"
#include "LDTest.h"
#include "MarkerRegression.h"
#include "HWETest.h"
#include "QTL.h"
#include "GeneInteractionBinary.h"
#include "GeneInteractionQTL.h"
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/erase.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#ifdef _WIN32
#include <windows.h>
#endif
#define VERSION 2.0

namespace po = boost::program_options;
std::stringstream report;

struct arguments {
  arguments()
      : haploAnalysis(false),
        assocAnalysis(false),
        html(true),
        hweAnalysis(false),
        ldAnalysis(false),
        epistasis(false),
        permutation(-1),
        webserver(false),
        qtl(false),
        lft(0.03),
        adjust(false),
        forceRegress(false),
        covar(""),
        epilb(2),
        epiub(2),
        epipermutation(50),
        model(SHEsis::ADDICTIVE),
        hapmethod(EM),
        HapLDT(0.1){};
  std::vector<std::string> inputfiles;
  std::vector<std::string> inputcases;
  std::vector<std::string> inputctrls;
  std::vector<std::string> snpnames;
  std::string output;
  std::string covar;
  int ploidy;
  int permutation;
  int epilb;
  int epiub;
  int epipermutation;
  bool containsPhenotype;
  bool haploAnalysis;
  bool assocAnalysis;
  bool qtl;
  bool hweAnalysis;
  bool ldAnalysis;
  bool webserver;
  bool adjust;
  bool forceRegress;
  bool epistasis;
  double lft;
  double HapLDT;
  std::vector<short> mask;
  SHEsis::LD_TYPE ldtype;
  HapMethod hapmethod;
  SHEsis::DiseaseModel model;
  bool html;
} SHEsisArgs;

void addOptions(int argc, char* argv[], po::options_description& desc,
                po::variables_map& vm);
void checkOptions(po::options_description& desc, po::variables_map& vm);
int ReadInput(int ploidy, bool containsPhenotype, std::string filepath,
              std::vector<std::vector<std::string> >& filecontent,std::vector<int>& invalidIdx);
boost::shared_ptr<SHEsis::SHEsisData> parseDataWithPhenotype(
    int snpnum, int ploidy, std::vector<std::vector<std::string> >& content);
boost::shared_ptr<SHEsis::SHEsisData> parseDataNoPhenotype(
    int snpnum, int ploidy, std::vector<std::vector<std::string> >& casecontent,
    std::vector<std::vector<std::string> >& ctrlcontent);
boost::shared_ptr<SHEsis::SHEsisData> parseInput();
std::vector< std::vector<double> > ReadCovar(std::string filepath,std::vector<int>& invalidIdx);
void getSnpNamefile(std::string& path, std::vector<std::string>& snp);
void getSnpNameline(std::string& names, std::vector<std::string>& snp);
void reportHtml(std::stringstream& report,
                boost::shared_ptr<SHEsis::AssociationTest> AssocHandle,
                boost::shared_ptr<SHEsis::MarkerRegression> RegressHandle,
                boost::shared_ptr<SHEsis::QTL> QTLHandle,
                boost::shared_ptr<SHEsis::HWETest> HWEHandle,
                boost::shared_ptr<SHEsis::HaplotypeLD> HapHandle,
                boost::shared_ptr<SHEsis::LDTest> LDHandle,
                boost::shared_ptr<SHEsis::GeneInteraction> EpiHandle);
void reporttxt(std::stringstream& report,
               boost::shared_ptr<SHEsis::AssociationTest> AssocHandle,
               boost::shared_ptr<SHEsis::MarkerRegression> RegressHandle,
               boost::shared_ptr<SHEsis::QTL> QTLHandle,
               boost::shared_ptr<SHEsis::HWETest> HWEHandle,
               boost::shared_ptr<SHEsis::HaplotypeLD> HapHandle,
               boost::shared_ptr<SHEsis::LDTest> LDHandle,
               boost::shared_ptr<SHEsis::GeneInteraction> EpiHandle);
void writePendingPage(std::ofstream& report);

int main(int argc, char* argv[]) {

  po::options_description desc_("Allowed options");
  po::variables_map vm;
  boost::shared_ptr<SHEsis::SHEsisData> data;
  addOptions(argc, argv, desc_, vm);
  try {
    checkOptions(desc_, vm);
    data = parseInput();
  }
  catch (std::runtime_error& e) {
    std::cout << "***ERROR: " << e.what() << "\n";
    std::cout << desc_ << "\n";
    exit(-1);
  };

  std::ofstream ofile;
  std::string filename =
      SHEsisArgs.output + (SHEsisArgs.html ? ".html" : ".txt");
  ofile.open(filename.c_str());
  if (!ofile) {
    std::cout << "***ERROR: cannot open output file " << filename << "\n";
    exit(-1);
  }
  if (SHEsisArgs.webserver) {
    writePendingPage(ofile);
  }
  ofile.close();

  if ((SHEsisArgs.snpnames.size() != 0) &&
      (data->getSnpNum() != SHEsisArgs.snpnames.size())) {
    std::cout << "***WARNING: given number of snp names is not the same at "
                 "that in data. Ignoring snp names...\n";
    SHEsisArgs.snpnames.clear();
  };
  for (int i = 0; i < SHEsisArgs.snpnames.size(); i++) {
    data->setLocusName(i, SHEsisArgs.snpnames[i]);
  }

  boost::shared_ptr<SHEsis::AssociationTest> AssocHandle;
  boost::shared_ptr<SHEsis::MarkerRegression> RegressHandle;
  boost::shared_ptr<SHEsis::QTL> QTLHandle;
  boost::shared_ptr<SHEsis::HWETest> HWEHandle;
  boost::shared_ptr<SHEsis::HaplotypeLD> HapHandle;
  boost::shared_ptr<SHEsis::LDTest> LDHandle;
  boost::shared_ptr<SHEsis::GeneInteraction> EpiHandle;

  if (SHEsisArgs.assocAnalysis) {
    std::cout << "Starting association test...\n";
    if(SHEsisArgs.forceRegress || SHEsisArgs.covar != ""){
    	RegressHandle.reset(new SHEsis::MarkerRegression(data));
    	if(SHEsisArgs.adjust)
    		RegressHandle->setAdjust(true);
        if (SHEsisArgs.permutation != -1)
        	RegressHandle->setPermutation(SHEsisArgs.permutation);
    	RegressHandle->regressAll();
    }else if(data->vQuantitativeTrait.size() == 0) {
      AssocHandle.reset(new SHEsis::AssociationTest(data));
      if (SHEsisArgs.adjust) AssocHandle->setAdjust(true);
      if (SHEsisArgs.permutation != -1) {
        AssocHandle->setPermutationTimes(SHEsisArgs.permutation);
        AssocHandle->permutation();
      } else {
        AssocHandle->association();
      }
    } else {
      QTLHandle.reset(new SHEsis::QTL(data));
      if (SHEsisArgs.adjust) QTLHandle->setAdjust(true);
      if (SHEsisArgs.permutation != -1) {
        QTLHandle->setPermutation(SHEsisArgs.permutation);
        QTLHandle->QTLPermutation();
      } else {
        QTLHandle->QTLTest();
      }
    }
    std::cout << "done\n";
  };

  if (SHEsisArgs.hweAnalysis && SHEsisArgs.ploidy > 1) {
    std::cout << "Starting Hardy-Weinberg equilibrium test...\n";
    HWEHandle.reset(new SHEsis::HWETest(data));
    HWEHandle->AllSnpHWETest();
    std::cout << "done\n";
  }
  if(SHEsisArgs.epistasis){
	  std::cout << "Starting epistasis analysis...\n";
	  if(data->vQuantitativeTrait.size() == 0){
		  EpiHandle.reset(new SHEsis::GeneInteractionBinary(data));
	  }else{
		  EpiHandle.reset(new SHEsis::GeneInteractionQTL(data));
//		  EpiHandle->setBinNum(SHEsisArgs.epibin);
//		  EpiHandle->setMinBin(SHEsisArgs.epiminbin);
//		  EpiHandle->setSamplePerBin(SHEsisArgs.epiminsample);
//		  std::cout<<"binNum="<<SHEsisArgs.epibin<<",minBin="<<SHEsisArgs.epiminbin<<",epiminisample="<<SHEsisArgs.epiminsample<<"\n";
	  }
	  EpiHandle->setAdjust(SHEsisArgs.adjust);
	  EpiHandle->setlb(SHEsisArgs.epilb);
	  EpiHandle->setub(SHEsisArgs.epiub);
	  EpiHandle->setPermutation(SHEsisArgs.epipermutation);
	  std::cout<<"epilb="<<SHEsisArgs.epilb<<",epiub="<<SHEsisArgs.epiub<<",epipermutation="<<SHEsisArgs.epipermutation<<"\n";
	  EpiHandle->CalGeneInteraction();
	  std::cout << "done\n";
  }

  if (SHEsisArgs.haploAnalysis && SHEsisArgs.ploidy > 1) {
    std::cout << "Starting haplotype analysis...\n";
    int snpnum = 0;
    if (SHEsisArgs.mask.size() != 0 &&
        SHEsisArgs.mask.size() != data->getSnpNum()) {
      std::cout << "***WARNING: length of given mask is wrong. Ignoring "
                   "mask...\n";
      SHEsisArgs.mask.clear();
    }
    for (int i = 0; i < SHEsisArgs.mask.size(); i++) {
      snpnum += SHEsisArgs.mask[i];
    };
    if (0 == snpnum)  // no mask
      HapHandle.reset(new SHEsis::HaplotypeLD(data));
    else
      HapHandle.reset(new SHEsis::HaplotypeLD(data, snpnum, SHEsisArgs.mask));

    if (SHEsisArgs.adjust) HapHandle->setAdjust(true);
    if (SHEsisArgs.lft >= 0 && SHEsisArgs.lft < 1)
      HapHandle->setFreqThreshold(SHEsisArgs.lft);
    else
      std::cout << "***WARNING: lowest frequency threshold for haplotype "
                   "analysis is invalid..defaulting to 0.03\n";
    HapHandle->setSilent(false);
    HapHandle->setLDT(SHEsisArgs.HapLDT);
    HapHandle->phaseAll();
    HapHandle->AssociationTest();
    std::cout << "done\n";
  }

  if (SHEsisArgs.ldAnalysis && SHEsisArgs.ploidy > 1) {
    std::cout << "Starting linkage disequilibrium analysis...\n";
    LDHandle.reset(new SHEsis::LDTest(data, SHEsisArgs.output));
    LDHandle->setLDType(SHEsisArgs.ldtype);
    LDHandle->AllLociLDtest();
    LDHandle->DrawLDMapDandR2();
    std::cout << "done\n";
  };


  if (SHEsisArgs.html)
    reportHtml(report, AssocHandle, RegressHandle,QTLHandle, HWEHandle, HapHandle, LDHandle,EpiHandle);
  else
    reporttxt(report, AssocHandle,RegressHandle, QTLHandle, HWEHandle, HapHandle, LDHandle,EpiHandle);
  ofile.open(filename.c_str());
  ofile << report.str();
  ofile.close();
  std::cout << "Results saved to " << filename << ".\n";
  return 0;
}

void writePendingPage(std::ofstream& report) {
  report << HtmlHeaderServer;
  report << "<head><meta http-equiv=\"refresh\" content=\"3\"></head>";
  report << "<div class=\"container\" style=\"padding:100px 100px 100px "
            "100px\">\n";
  report << "<div style=\"line-height:40px\">\n";
  report << "<font size=\"5\">\n";
  report << "<p>We have received your request. This page will auto-fresh until "
            "the results are ready.</p>";
  report << "<p>The notification of completed result will also be sent to the "
            "E-mail address you provided.</p>";
  report << "<p>The results will be available for 72 hours in web, after which "
            "they will be deleted.</p>";
  report << "</font></div></div></div></body></html>\n";
}

void reporttxt(std::stringstream& report,
               boost::shared_ptr<SHEsis::AssociationTest> AssocHandle,
               boost::shared_ptr<SHEsis::MarkerRegression> RegressHandle,
               boost::shared_ptr<SHEsis::QTL> QTLHandle,
               boost::shared_ptr<SHEsis::HWETest> HWEHandle,
               boost::shared_ptr<SHEsis::HaplotypeLD> HapHandle,
               boost::shared_ptr<SHEsis::LDTest> LDHandle,
               boost::shared_ptr<SHEsis::GeneInteraction> EpiHandle) {
  if (AssocHandle) {
    report << AssocHandle->reporttxt();
  }
  if (QTLHandle) {
    report << QTLHandle->reporttxt();
  }
  if(RegressHandle){
	  report << RegressHandle->reporttxt();
  }
  if (HWEHandle) {
    report << HWEHandle->reporttxt();
  }
  if (HapHandle) {
    report << HapHandle->reporttxt();
  }
  if (LDHandle) {
    report << LDHandle->reporttxt();
  }
  if(EpiHandle){
	  report<<EpiHandle->reporttxt();
  }
}

void reportHtml(std::stringstream& report,
                boost::shared_ptr<SHEsis::AssociationTest> AssocHandle,
                boost::shared_ptr<SHEsis::MarkerRegression> RegressHandle,
                boost::shared_ptr<SHEsis::QTL> QTLHandle,
                boost::shared_ptr<SHEsis::HWETest> HWEHandle,
                boost::shared_ptr<SHEsis::HaplotypeLD> HapHandle,
                boost::shared_ptr<SHEsis::LDTest> LDHandle,
                boost::shared_ptr<SHEsis::GeneInteraction> EpiHandle) {
  if (!SHEsisArgs.webserver) {
    report << HtmlHeader;
    report << "<h1>SHEsis </h1>\n";
  } else {
    report << HtmlHeaderServer;
  }
  if (AssocHandle) {
    report << AssocHandle->reporthtml();
  }
  if (QTLHandle) {
    report << QTLHandle->reporthtml();
  }
  if( RegressHandle){
	  report<< RegressHandle->reporthtml();
  }
  if (HWEHandle) {
    report << HWEHandle->reporthtml();
  }
  if(EpiHandle){
	  report<<EpiHandle->reporthtml();
  }
  if (HapHandle) {
    report << HapHandle->reporthtml();
  }
  if (LDHandle) {
    report << LDHandle->reporthtml();
  }

  if (!SHEsisArgs.webserver)
    report << "</body>\n</html>\n";
  else
    report << "</div></body>\n</html>\n";
}

boost::shared_ptr<SHEsis::SHEsisData> parseInput() {
  boost::shared_ptr<SHEsis::SHEsisData> data;
  std::vector<std::vector<std::string> > filecontent;
  std::vector<int> invalidIdx;
  std::vector<std::vector<std::string> > filecontentcase;
  std::vector<std::vector<std::string> > filecontentctrl;
  int snpnum, snpcase, snpctrl;
  if (SHEsisArgs.inputfiles.size() > 0) {
	BOOST_ASSERT(SHEsisArgs.inputfiles.size()==1);
    for (int i = 0; i < SHEsisArgs.inputfiles.size(); i++) {
      snpnum = ReadInput(SHEsisArgs.ploidy, SHEsisArgs.containsPhenotype,
                         SHEsisArgs.inputfiles[i], filecontent,invalidIdx);
    }
    data = parseDataWithPhenotype(snpnum, SHEsisArgs.ploidy, filecontent);
  } else if (SHEsisArgs.inputctrls.size() > 0 &&
             SHEsisArgs.inputcases.size() > 0) {
	BOOST_ASSERT(SHEsisArgs.inputctrls.size()==1 && SHEsisArgs.inputcases.size()==1 );
    for (int i = 0; i < SHEsisArgs.inputctrls.size(); i++)
      snpctrl = ReadInput(SHEsisArgs.ploidy, SHEsisArgs.containsPhenotype,
                          SHEsisArgs.inputctrls[i], filecontentctrl,invalidIdx);
    for (int i = 0; i < SHEsisArgs.inputcases.size(); i++)
      snpcase = ReadInput(SHEsisArgs.ploidy, SHEsisArgs.containsPhenotype,
                          SHEsisArgs.inputcases[i], filecontentcase,invalidIdx);
    if (snpctrl != snpcase)
      throw std::runtime_error("SNP number in cases and controls disagrees.");
    data = parseDataNoPhenotype(snpcase, SHEsisArgs.ploidy, filecontentcase,
                                filecontentctrl);
  } else {
    throw std::runtime_error("Please check the input files.");
  };
  if(SHEsisArgs.covar!=""){
	  data->covar=ReadCovar(SHEsisArgs.covar,invalidIdx);
  }

  return data;
}

boost::shared_ptr<SHEsis::SHEsisData> parseDataWithPhenotype(
    int snpnum, int ploidy, std::vector<std::vector<std::string> >& content) {
  int samplenum = content.size();
  boost::shared_ptr<SHEsis::SHEsisData> pdata(
      new SHEsis::SHEsisData(samplenum, snpnum, ploidy));
  for (int sample = 0; sample < samplenum; sample++) {
	  pdata->sampleName.push_back(content[sample][0].c_str());
    if (SHEsisArgs.qtl) {
      pdata->vQuantitativeTrait.push_back(atof(content[sample][1].c_str()));
    } else {
      pdata->vLabel.push_back(
          (SHEsis::SampleStatus)std::atoi(content[sample][1].c_str()));
    }
    for (int snp = 0; snp < snpnum; snp++) {
      for (int p = 0; p < ploidy; p++) {
        pdata->mGenotype[sample][snp][p] =
            pdata->GetAlleleCode(content[sample][snp * ploidy + 2 + p]);
      }
    }
  };
  return pdata;
}

boost::shared_ptr<SHEsis::SHEsisData> parseDataNoPhenotype(
    int snpnum, int ploidy, std::vector<std::vector<std::string> >& casecontent,
    std::vector<std::vector<std::string> >& ctrlcontent) {
  int casenum = casecontent.size();
  int ctrlnum = ctrlcontent.size();
  boost::shared_ptr<SHEsis::SHEsisData> pdata(
      new SHEsis::SHEsisData(casenum + ctrlnum, snpnum, ploidy));
  for (int sample = 0; sample < casecontent.size(); sample++) {
    pdata->vLabel.push_back(SHEsis::CASE);
    for (int snp = 0; snp < snpnum; snp++) {
      for (int p = 0; p < ploidy; p++) {
        pdata->mGenotype[sample][snp][p] =
            pdata->GetAlleleCode(casecontent[sample][snp * ploidy + 1 + p]);
      }
    }
  };
  for (int sample = 0; sample < ctrlcontent.size(); sample++) {
    pdata->vLabel.push_back(SHEsis::CONTROL);
    for (int snp = 0; snp < snpnum; snp++) {
      for (int p = 0; p < ploidy; p++) {
        pdata->mGenotype[sample + casenum][snp][p] =
            pdata->GetAlleleCode(ctrlcontent[sample][snp * ploidy + 1 + p]);
      }
    }
  };
  return pdata;
}

void addOptions(int argc, char* argv[], po::options_description& desc,
                po::variables_map& vm) {
  desc.add_options()("help", "produce help message")(
      "input", po::value<std::vector<std::string> >(),
      "path for the input file containing both cases and controls"
      /*", can be specified for multiple times"*/)(
      "input-case", po::value<std::vector<std::string> >(),
      "path for the input file containing cases "
      /*", can be specified for multiple times"*/)("input-ctrl", po::value<std::vector<std::string> >(),
               "path for the input file containing controls "
               /*", can be specified for multiple times"*/)(
      "covar",po::value<std::string>(),"path for covariates")(
      "snpname-file", po::value<std::string>(),
      "path for file that contains names of snps")(
      "snpname-line", po::value<std::string>(), "snp names are as arguments")(
      "output", po::value<std::string>(), "prefix of output files")(
      "report-txt",
      "report results in plain-text format. By default, results will be "
      "reported in html.")("ploidy", po::value<int>(), "number of ploidy")(
      "hwe", "perform Hardy-Weinberg disequilibrium test")(
      "assoc",
      "perform association test, case/control analysis by default. To perform "
      "quantitative trait loci analysis, please specified together with "
      "--qtl.")("qtl",
                "input phenotype is quantitative traits. input file should be "
                "specified with --input, the second column of the input file "
                "is the quantitative trait")("permutation", po::value<int>(),
                                             "times for permutation")(
      "regress","force regression, even if no covariates provided.")(
       "model",po::value<std::string>(),"diseas model for regression analysis. possible value:dominant|reccesive|addictive,"
    		   "default:addictive")(
      "haplo-EM",
      "perform haplotype analysis using expectation maximization algorithm")(
      "haplo-SAT", "perform haplotype analysis using SAT-based algorithm, no longer used")(
      "mask", po::value<std::string>(),
      "mask of snps for haplotype analysis, comma delimited. eg. mask=1,0,1 to "
      "use 1st and 3rd SNPs when there are 3 SNPs in all.")(
      "lft", po::value<double>(),
      "lowest frequency threshold for haplotype analysis")(
      "hapld",po::value<double>(),
      "R2 to define a LD block. Snp in one snp block will be phased together. default:0.2")(
      "epistasis","perform gene interaction analysis")(
      "epi-lb",po::value<int>(),"lower bound of number of snps to perform gene interaction analysis,default:2 (2-way interaction)")(
      "epi-ub",po::value<int>(),"upper bound of number of snps to perform gene interaction analysis,default:2 (2-way interaction)")(
      "epi-permutation",po::value<int>(),"number of permutations when performing gene interaction analysis, default:50")(
      "ld-in-case", "perform Linkage disequilibrium test in cases")(
      "ld-in-ctrl", "perform Linkage disequilibrium test in controls")(
      "ld", "perform Linkage disequilibrium test in both cases and controls")(
      "adjust", "adjust p-value for multiple testing")(
      "webserver", "Internal use for webserver");
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
}

void getSnpNamefile(std::string path, std::vector<std::string>& snp) {
  std::string line;
  std::vector<std::string> strs;
  std::ifstream file(path.c_str());
  if (!file.is_open()) throw std::runtime_error("Cannot open snp name file.");
  while (getline(file, line)) {
    if (line.empty()) continue;
    boost::trim_if(line, boost::is_any_of("\t "));
    boost::erase_all(line, "\r");
    boost::split(strs, line, boost::is_any_of("\t "), boost::token_compress_on);
    if (strs.size() > 1) {
      std::cout << "***WARNING: " << path
                << " contains more than 1 fileds in a single line.";
      std::cout << "Ignoring this file...\n";
      snp.clear();
      return;
    }
    snp.push_back(line);
  }
}

void getSnpNameline(std::string names, std::vector<std::string>& snp) {
  if (names.empty()) return;
  boost::trim_if(names, boost::is_any_of("\t ,"));
  boost::split(snp, names, boost::is_any_of("\t ,"), boost::token_compress_on);
}

void checkOptions(po::options_description& desc, po::variables_map& vm) {
  if (vm.count("help")) {
    std::cout << desc << "\n";
    exit(0);
  }

  if ((0 == vm.count("input")) && (0 == vm.count("input-case")) &&
      (0 == vm.count("input-ctrl")))
    throw std::runtime_error("no input file specified.");
  if (vm.count("input") && (vm.count("input-case") || vm.count("input-ctrl")))
    throw std::runtime_error(
        "--input and --input-case/--input-ctrl cannot be specified at the same "
        "time.");
  if (vm.count("qtl") > 0 && (vm.count("input-case") || vm.count("input-ctrl")))
    throw std::runtime_error(
        "input file for quantitative trait analaysis should be specified with "
        "--input.");
  if (vm.count("qtl") > 0 && (vm.count("input") == 0))
    throw std::runtime_error("no input file specified.");
  if ((0 == vm.count("input-case") && (0 != vm.count("input-ctrl"))) ||
      (0 == vm.count("input-ctrl") && (0 != vm.count("input-case"))))
    throw std::runtime_error(
        "--input-case and --input-ctrl should both be specified.");
  if (vm.count("input") > 0) {
    SHEsisArgs.inputfiles = vm["input"].as<std::vector<std::string> >();
    SHEsisArgs.containsPhenotype = true;
  } else if (vm.count("input-case") > 0 && vm.count("input-ctrl") > 0) {
    SHEsisArgs.inputcases = vm["input-case"].as<std::vector<std::string> >();
    SHEsisArgs.inputctrls = vm["input-ctrl"].as<std::vector<std::string> >();
    SHEsisArgs.containsPhenotype = false;
  } else {
    throw std::runtime_error("error in parsing input files.");
  }
  if (0 != vm.count("snpname-file") && 0 != vm.count("snpname-line"))
    throw std::runtime_error(
        "--snpname-line and --snpname-file cannot be specified at the same "
        "time.");
  if (vm.count("snpname-file")) {
    getSnpNamefile(vm["snpname-file"].as<std::string>(), SHEsisArgs.snpnames);
  } else if (vm.count("snpname-line")) {
    getSnpNameline(vm["snpname-line"].as<std::string>(), SHEsisArgs.snpnames);
  }
  if( 0 != vm.count("covar"))
	SHEsisArgs.covar=vm["covar"].as<std::string>();
  if (0 == vm.count("output"))
    SHEsisArgs.output = "SHEsis";
  else
    SHEsisArgs.output = vm["output"].as<std::string>();

  if (0 != vm.count("lft")) SHEsisArgs.lft = vm["lft"].as<double>();

  if (0 == vm.count("ploidy"))
    throw std::runtime_error("no ploidy information given.");
  else if (vm["ploidy"].as<int>() < 1)
    throw std::runtime_error("number of ploidy should be higher than 1.");
  else
    SHEsisArgs.ploidy = vm["ploidy"].as<int>();

  int ldcase = vm.count("ld-in-case");
  int ldctrl = vm.count("ld-in-ctrl");
  int ld = vm.count("ld");
  if (ldcase + ldctrl + ld > 1) {
    throw std::runtime_error(
        "--ld-in-case/--ld-in-ctrl/--ld cannot be specified at the same time");
  } else if (ldcase + ldctrl + ld == 1) {
    SHEsisArgs.ldAnalysis = true;
    if (ldcase) SHEsisArgs.ldtype = SHEsis::LD_IN_CASE;
    if (ldctrl) SHEsisArgs.ldtype = SHEsis::LD_IN_CTRL;
    if (ld) SHEsisArgs.ldtype = SHEsis::LD_IN_BOTH;
  };

//  if (vm.count("permutation") > 0 && vm.count("assoc") == 0)
  //  throw std::runtime_error("--permutaion should be used along with --assoc");
  //  if(vm.count("assoc")!=0 && vm.count("qtl")!=0)
  //	  throw std::runtime_error("--assoc cannot be specified together with
  //--qtl");
  if (vm.count("assoc")) {
    SHEsisArgs.assocAnalysis = true;
    if (vm.count("permutation")) {
      SHEsisArgs.permutation = vm["permutation"].as<int>();
      if (SHEsisArgs.permutation <= 0)
        throw std::runtime_error("permutation times should be higher than 0.");
    }
  };
  if (vm.count("qtl")) {
    SHEsisArgs.qtl = true;
  }
  if (vm.count("model")) {
    std::string m = vm["model"].as<std::string>();
    if(m == "addictive")
    	SHEsisArgs.model=SHEsis::ADDICTIVE;
    else if (m == "dominant")
    	SHEsisArgs.model=SHEsis::DOMINANT;
    else if (m== "recessive")
    	SHEsisArgs.model=SHEsis::RECESSIVE;
    else
    	throw std::runtime_error("disease model should be addictive|dominant|recessive.");
  }
  if (vm.count("regress")) {
    SHEsisArgs.forceRegress = true;
  }

  if (vm.count("hwe")) {
    SHEsisArgs.hweAnalysis = true;
  };
  if (vm.count("haplo-SAT") != 0 && vm.count("haplo-EM") != 0) {
    throw std::runtime_error(
        "--haplo-SAT and --haplo-EM cannot be specified at the same time.");
  };
  if (vm.count("haplo-SAT") != 0) {
    SHEsisArgs.haploAnalysis = true;
    SHEsisArgs.hapmethod = SAT;
  };
  if (vm.count("haplo-EM") != 0) {
    SHEsisArgs.haploAnalysis = true;
    SHEsisArgs.hapmethod = EM;
  };
  if( vm.count("hapld")!=0){
	  SHEsisArgs.HapLDT=vm["hapld"].as<double>();
  }

  if(vm.count("epistasis") !=0){
	  SHEsisArgs.epistasis=true;
  }
  if(vm.count("epi-lb")!=0){
	  SHEsisArgs.epilb = vm["epi-lb"].as<int>();
	  if(SHEsisArgs.epilb<2)
		  throw std::runtime_error("--epi-lb should be above 2");
  }
  if(vm.count("epi-ub")!=0){
	  SHEsisArgs.epiub = vm["epi-ub"].as<int>();
	  if(SHEsisArgs.epiub<SHEsisArgs.epilb)
		  throw std::runtime_error("--epi-ub should be equal or higher than --epi-lb");
  }
  if(vm.count("epi-permutation")!=0){
	  SHEsisArgs.epipermutation= vm["epi-permutation"].as<int>();
	  if(SHEsisArgs.epipermutation<=0)
		  throw std::runtime_error("--epi-permutation should be higher than 0");
  }

  if (vm.count("adjust") != 0) SHEsisArgs.adjust = true;

  if (vm.count("mask")) {
    std::string maskstr = vm["mask"].as<std::string>();
    std::vector<std::string> maskvec;
    boost::split(maskvec, maskstr, boost::is_any_of("\t ,"));
    for (int i = 0; i < maskvec.size(); i++) {
      if (std::strcmp("0", maskvec[i].c_str()) == 0) {
        SHEsisArgs.mask.push_back(0);
      } else if (std::strcmp("1", maskvec[i].c_str()) == 0) {
        SHEsisArgs.mask.push_back(1);
      } else {
        // throw std::runtime_error("mask should be composed of 0 and 1. But " +
        //                        maskvec[i] + " found.");
        SHEsisArgs.mask.clear();
        std::cout << "***WARNING: mask should be composed of 0 and 1. But " +
                         maskvec[i] + " found. No mask will be used. \n";
      }
    }
  }
  if (vm.count("webserver")) {
    SHEsisArgs.webserver = true;
  }
  if (vm.count("report-txt")) SHEsisArgs.html = false;

  if (!SHEsisArgs.hweAnalysis && !SHEsisArgs.haploAnalysis &&
      !SHEsisArgs.assocAnalysis && !SHEsisArgs.ldAnalysis && !SHEsisArgs.epistasis)
    throw std::runtime_error(
        "at least one type of analysis should be specified.");
}

std::vector< std::vector<double> > ReadCovar(std::string filepath,std::vector<int>& invalidIdx){
	std::vector< std::vector<double> > res;
	std::string line;
	std::ifstream file(filepath.c_str());
	if (!file.is_open())
	   throw std::runtime_error("Cannot open input file: " + filepath);
	int lineidx=1;
	while (getline(file,line)){
		std::vector<std::string> strs;
		boost::erase_all(line,"\r");
		boost::trim_if(line,boost::is_any_of("\t ,"));
		boost::split(strs,line,boost::is_any_of("\t ,"),boost::token_compress_on);
		std::vector<double> aCovar;
		for(int i=0;i<strs.size();i++){
			double c;
			try{
				c=boost::lexical_cast<double>(strs[i]);
			}catch(const boost::bad_lexical_cast &){
				if(strs[i]=="NA"){
					c=DBL_MAX;
				}else{
					std::stringstream ss;
					ss<<"Error in line "<<lineidx<<", file:"<<filepath<<
							",covariates should be numerical. For missing data, please use 'NA'. Found "<<strs[i]<<".";
					throw std::runtime_error(ss.str());
				};

			}
			aCovar.push_back(c);
		}
		if(res.size()!=0 && res[0].size()!=aCovar.size()){
			std::stringstream ss;
			ss<<"Error in line "<<lineidx<<", file:"<<filepath<<
					","<<res[0].size()<<" fileds expected, but found "<<aCovar.size()<<".";
			throw std::runtime_error(ss.str());
		}
		res.push_back(aCovar);
	}
	for(int i=0;i<invalidIdx.size();i++){
		res.erase(res.begin()+invalidIdx[i]);
	}
	return res;
}

int ReadInput(int ploidy, bool containsPhenotype, std::string filepath,
              std::vector<std::vector<std::string> >& filecontent,std::vector<int>& invalidIdx) {
  std::string line;
  std::ifstream file(filepath.c_str());
  if (!file.is_open())
    throw std::runtime_error("Cannot open input file: " + filepath);
  int lineidx = 1;
  int snpnum = 0;
  int expectedfileds = 0;
  while (getline(file, line)) {
    std::vector<std::string> strs;
    boost::erase_all(line, "\r");
    boost::trim_if(line, boost::is_any_of("\t ,"));
    if (line.empty()) {
      lineidx++;
      continue;
    }
    boost::split(strs, line, boost::is_any_of("\t ,"),
                 boost::token_compress_on);

    if (containsPhenotype) {
    	if("NA" == strs[1]){
    		invalidIdx.push_back(filecontent.size());
    		continue;
    	}
    	if(!SHEsisArgs.qtl){
    		if(strs[1] != "case"&& strs[1] != "ctrl"){
    			int pheno;
				try{
					pheno=boost::lexical_cast<int>(strs[1]);
				}catch(const boost::bad_lexical_cast &){
					std::stringstream ss;
					ss<<"Error in line "<<lineidx<<", file:"<<filepath<<",phenotype should be 1 or ctrl for controls, 2 or case for cases, but found "<<strs[1]<<".";
					throw std::runtime_error(ss.str());
				}
				if(pheno!=1 && pheno!=2){
					std::stringstream ss;
					ss<<"Error in line "<<lineidx<<", file:"<<filepath<<",phenotype should be 1 or ctrl for controls, 2 or case for cases, but found "<<strs[1]<<".";
					throw std::runtime_error(ss.str());
				}
    		}else{

    			strs[1]=(strs[1]=="case"?"2":"1");
    		}
    	}else{
			try{
				int pheno=boost::lexical_cast<float>(strs[1]);
			}catch(const boost::bad_lexical_cast &){
				std::stringstream ss;
				ss<<"Error in line "<<lineidx<<", file:"<<filepath<<",phenotype should be numerical, but found "<<strs[1]<<".";
				throw std::runtime_error(ss.str());
			}
    	}
      if (!SHEsisArgs.qtl &&
          (SHEsis::SampleStatus)(std::atoi(strs[1].c_str()) != SHEsis::CASE &&
                                 (SHEsis::SampleStatus)(std::atoi(
                                     strs[1].c_str())) != SHEsis::CONTROL)) {
        std::stringstream ss;
        ss << "Error in line " << lineidx << ", file:" << filepath
           << ", phenotype should be either " << SHEsis::CASE
           << " for cases or " << SHEsis::CONTROL << " for controls, but "
           << strs[1]
           << " found. For quantitative trait analysis, please use --qtl.";
        throw std::runtime_error(ss.str());
      };
    };

    if (expectedfileds == 0) {
      if ((strs.size() - 1 - (int)containsPhenotype) % ploidy != 0)
        throw std::runtime_error("Error in line 1, file:" + filepath +
                                 ". Either snp num or ploidy num is wrong");
      snpnum = (strs.size() - 1 - (int)containsPhenotype) / ploidy;
      expectedfileds = (1 + (int)containsPhenotype + snpnum * ploidy);
      if (filecontent.size() != 0 && filecontent[0].size() != expectedfileds) {
        throw std::runtime_error("number of snps disagree between files");
      }
    }

    if (strs.size() != expectedfileds) {
      std::stringstream ss;
      ss << "Error in line " << lineidx << ", file:" << filepath
         << ", expecting " << expectedfileds << " fileds, but " << strs.size()
         << " fileds are found";
      throw std::runtime_error(ss.str());
    };
    lineidx++;
    filecontent.push_back(strs);
    if (filecontent.size() == 0)
      throw std::runtime_error("Empty file: " + filepath);
  }
  if (snpnum == 0) throw std::runtime_error("No snp data found in " + filepath);
  return snpnum;
}
