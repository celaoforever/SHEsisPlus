/*
 * Haplotype.cpp
 *
 *  Created on: Aug 19, 2014
 *      Author: ada
 */
#include "fisher.h"
#include "utility.h"
#include "Haplotype.h"
#include <math.h>
#include <sstream>
#include "System.h"
#include "ParseUtils.h"
#include "Options.h"
#include "Dimacs.h"
#include "Solver.h"
#include <boost/lexical_cast.hpp>


namespace SHEsis {
#define CEIL(x) (x-(int)x)>0?((int)x+1):(x)
#define MIN(a,b) a<=b?a:b
#define ABS(a) a<0?(-1*a):a
#define EOL " 0\n";this->ClauseNum++;
#define EOL2 "0\n";this->ClauseNum++;

void Haplotype::statOccurence(){
	for(int iSample=0;iSample<data->getSampleNum();iSample++){
		for(int iSnp=0;iSnp<data->getSnpNum();iSnp++){
			this->occurence[iSample][iSnp].resize(data->vLocusInfo[iSnp].BothAlleleCount.size(),0);
			this->missing[iSample][iSnp]=0;
			for(int p=0;p<data->getNumOfChrSet();p++){
				if(GENOTYPE_MISSING != data->mGenotype[iSample][iSnp][p]){
					int idx=data->vLocusInfo[iSnp].getAlleleIndex(data->mGenotype[iSample][iSnp][p]);
					BOOST_ASSERT(-1 != idx);
					this->occurence[iSample][iSnp][idx]++;
				}else
				{
					this->missing[iSample][iSnp]++;
				}
			}
		}
	}
//	std::cout<<"\n";
//	for(int iSample=0;iSample<data->getSampleNum();iSample++){
//	for(int iSnp=0;iSnp<this->occurence.shape()[1];iSnp++){
//	for(int k=00;k<this->occurence[iSample][iSnp].size();k++){
//	std::cout<<this->occurence[iSample][iSnp][k]<<"/";
//	}
//	std::cout<<" ";
//	}
//	std::cout<<"\n";
//	};
}


void Haplotype::statOccurenceMask(){
//	std::cout<<"mask:";
//	for(int i=0;i<mask.size();i++){
//		std::cout<<mask[i]<<",";;
//	}
	std::cout<<"\n";
	int subiSnp;
	for(int iSample=0;iSample<data->getSampleNum();iSample++){
		subiSnp=0;
		for(int iSnp=0;iSnp<data->getSnpNum();iSnp++){
			if(mask[iSnp]){
				if(this->SnpIdx.size()==0 || iSnp>this->SnpIdx[SnpIdx.size()-1])
					this->SnpIdx.push_back(iSnp);
				this->occurence[iSample][subiSnp].resize(data->vLocusInfo[iSnp].BothAlleleCount.size(),0);
				this->missing[iSample][subiSnp]=0;
				for(int p=0;p<data->getNumOfChrSet();p++){
					if(GENOTYPE_MISSING != data->mGenotype[iSample][iSnp][p]){
						int idx=data->vLocusInfo[iSnp].getAlleleIndex(data->mGenotype[iSample][iSnp][p]);
						BOOST_ASSERT(-1 != idx);
						this->occurence[iSample][subiSnp][idx]++;
					}else{
						this->missing[iSample][subiSnp]++;
					}
				}
				subiSnp++;
			}
		}
	}
//	std::cout<<"\n";
//	for(int i=0;i<this->SnpIdx.size();i++){
//		std::cout<<SnpIdx[i]<<",";
//	}
//	std::cout<<"\n";
//	for(int iSample=0;iSample<data->getSampleNum();iSample++){
//	for(int iSnp=0;iSnp<this->occurence.shape()[1];iSnp++){
//	for(int k=00;k<this->occurence[iSample][iSnp].size();k++){
//	std::cout<<this->occurence[iSample][iSnp][k]<<"/";
//	}
//	std::cout<<" ";
//	}
//	std::cout<<"\n";
//	};


}
int Haplotype::solve(){
	//std::stringstream sat;

	Minisat::Solver S;
	parse_DIMACS(this->res, S);
    if (!S.simplify()){
        return 0;
    }
    Minisat::vec<Minisat::Lit> dummy;
    Minisat::lbool ret = S.solveLimited(dummy);
    if (ret == l_True){
//            sat<<"SAT\n";
            for (int i = 0; i < S.nVars(); i++)
                if (S.model[i] != l_Undef){
                	sat=sat+((i==0)?"":" ")+((S.model[i]==l_True)?"":"-")+boost::lexical_cast<std::string>(i+1);
//                    printf("%s%s%d", (i==0)?"":" ", (S.model[i]==l_True)?"":"-", i+1);}
                };
//           sat<<" 0\n";
//           std::cout<<sat;
           return 1;
        };
    return 0;

}

void Haplotype::parseSolution(IndexingVariables& variables,int assumed_haplotypes){
    std::vector<std::string> strs;
    std::stringstream tmpss;
    boost::split(strs,this->sat,boost::is_any_of("\t "));
    boost::shared_ptr<int[]> parities(new int[variables.numberOfVariables()+1]);
    parities[0]=variables.numberOfVariables();
    for(int i=1;i<=parities[0];i++){
    	parities[i]=0;
    };
//    std::cout<<"\n\n\nstrs:\n";
//    for(int i=0;i<strs.size();i++){
//    	std::cout<<strs[i].c_str()<<" ";
//    }
//	std::cout<<"\n\n\n";
    for(int i=0;i<strs.size();i++){
    	int var=boost::lexical_cast<int>(strs[i].c_str());
    	BOOST_ASSERT(0!=var);
    	int index=variables.enumeration[ABS(var)];
//    	std::cout<<index<<"\n";
    	BOOST_ASSERT(index<=parities[0]);
    	parities[index-1]=var;
    }
    variables.setParities(parities);
//
    for(int i=0;i<assumed_haplotypes;i++){
    	boost::shared_ptr<short[]> haplo(new short[this->SnpIdx.size()]);
    	for(int j=0;j<this->SnpIdx.size();j++){
    		int width=CEIL(log2(this->occurence[0][j].size()));
    		char* chr_=new char[width+1];
    		char* chr =new char[width+1];
    		chr_[width]='\0';
    		chr[width]='\0';
    		for(int k=0;k<width;k++){
    			tmpss.str("");
    			tmpss<<j<<"_"<<k<<"haplotypes";
    			int f=variables.getEvalutatedId(tmpss.str(),SetSharedPtr(1,i));
    			BOOST_ASSERT(0 != f);
    			if(f<0)
    				chr_[k]='0';
    			else
    				chr_[k]='1';
    		};
    		int counter=width-1;
    		for(int k=0;k<width;k++){
    			chr[counter--]=chr_[k];
    		}
    		std::string str_chr(chr);
    		int index=std::atoi(str_chr.c_str());
    		haplo[j]=this->data->vLocusInfo[this->SnpIdx[j]].getAlleleType(index);
    		delete[] chr;
    		delete[] chr_;
    	}
    	this->Results.haplotypes.push_back(haplo);
    	for(int i=0;i<this->SnpIdx.size();i++){
        	std::cout<<haplo[i];
    	};
    	std::cout<<"\n";

    }

    for(int i=0;i<this->data->getSampleNum();i++){
    	for(int k=0;k<this->data->getNumOfChrSet();k++){
    		this->Results.genotypes[i][k]=-1;
    		for(int j=0;j<assumed_haplotypes;j++){
    			int f=variables.getEvalutatedId("selections",SetSharedPtr(3,k,j,i));
    			if(f>0){
    				this->Results.genotypes[i][k]=j;
    				break;
    			}
    		}
    	}
    }

    for(int i=0;i<this->data->getSampleNum();i++){
    	for(int k=0;k<this->data->getNumOfChrSet();k++){
    		std::cout<<Results.genotypes[i][k]<<",";
    	}
    	std::cout<<"\n";
    };


}

void Haplotype::associationTest(){
	int haploNum=this->Results.haplotypes.size();
	this->Results.CaseCount.reset(new int[haploNum]);
	this->Results.ControlCount.reset(new int[haploNum]);
	for(int i=0;i<this->Results.haplotypes.size();i++){
		this->Results.CaseCount[i]=0;
		this->Results.ControlCount[i]=0;
	}
    for(int i=0;i<this->data->getSampleNum();i++){
    	for(int k=0;k<this->data->getNumOfChrSet();k++){
			int idx=Results.genotypes[i][k];
    		if(CASE == this->data->vLabel[i]){
    			this->Results.CaseCount[idx]++;
    		}else if(CONTROL == this->data->vLabel[i]){
    			this->Results.ControlCount[idx]++;
    		}
    	}
    }
    std::cout<<"contigency:\n";
    double* contigency= new double[2*haploNum];
    int idx=0;
    for(int i=0;i<haploNum;i++){
    	contigency[idx++]=this->Results.ControlCount[i];
    	contigency[idx++]=this->Results.CaseCount[i];
    	std::cout<<this->Results.ControlCount[i]<<","<<this->Results.CaseCount[i]<<"\n";
    };


	  double expect = -1.0;
	  double percnt = 100.0;
	  double emin = 0;
	  double pre = 0, prt = 0;
	  int ws = 300000;
	  int nrow=2;
	  fexact(&nrow, &haploNum, contigency, &nrow, &expect, &percnt, &emin, &prt, &pre, &ws);
	  this->Results.FisherP=pre;

	  //Pearson's ChiSquare test
	  PearsonChiSquareTest(contigency,nrow,haploNum,this->Results.ChiSquare,this->Results.PearsonP);
	  std::cout<<"fisherp:"<<this->Results.FisherP;
	  std::cout<<"\npearsonp:"<<this->Results.PearsonP;




}

void Haplotype::BuildModel(IndexingVariables& variables, int number_of_explaining_haplotypes){
//#define res std::cout
	std::stringstream tmpss;
	boost::multi_array<int,1> anti_haplotypes;
	int number_of_genotypes=this->occurence.shape()[0];
	int length_of_genotypes=this->occurence.shape()[1];
	int ploidy=this->data->getNumOfChrSet();
	int number_of_known_haplotypes=0;

	//IndexingVariables variables;
	createVariables(number_of_explaining_haplotypes,variables);
	for(int which_genotype=0;which_genotype<number_of_genotypes;which_genotype++){
		for(int which_index=0;which_index<length_of_genotypes;which_index++){
			int width=CEIL(log2(this->occurence[which_genotype][which_index].size()));

			if(0 == this->missing[which_genotype][which_index]){
				bool b=true;
				for(int i=0;i<this->occurence[which_genotype][which_index].size();i++){
					if(ploidy != this->occurence[which_genotype][which_index][i])
						continue;
					boost::shared_ptr<int[] >number=toBooleanInt(width,i);
					for(int which_chromosome=0;which_chromosome<ploidy;which_chromosome++){
						for(int k=0;k<width;k++){
							int f1=(int)pow(-1,1-number[k]);
							tmpss.str("");
							tmpss<<which_index<<"_"<<k<<"combinations";
							int f2=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,which_genotype));
							res<<(f1*f2)<<EOL;
						}
					}
					b=false;
					break;
				}
				if(b){
					getGeneralCoding(ploidy,which_genotype,which_index,variables);
				}
			}else if(this->missing[which_genotype][which_index]<ploidy){
					getGeneralCodingMissing(ploidy,which_genotype,which_index,variables);
			}else{
					getGeneralCodingTotalyMissing(ploidy,which_genotype,which_index,variables);
			}

			for(int which_explaining_genotype=0;which_explaining_genotype<number_of_explaining_haplotypes;which_explaining_genotype++){
				for(int which_chromosome=0;which_chromosome<ploidy;which_chromosome++){
					for(int k=0;k<width;k++){
						tmpss.str("");
						tmpss<<which_index<<"_"<<k<<"haplotypes";
						int f1=variables.getEnumeration(tmpss.str(),SetSharedPtr(1,which_explaining_genotype));
						tmpss.str("");
						tmpss<<which_index<<"_"<<k<<"combinations";
						int f2=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,which_genotype));
						int f3=variables.getEnumeration("selections",SetSharedPtr(3,which_chromosome,which_explaining_genotype,which_genotype));
						res<<f1<<" "<<(-1*f2)<<" "<<(-1*f3)<<EOL;
						tmpss.str("");
						tmpss<<which_index<<"_"<<k<<"haplotypes";
						int f4=variables.getEnumeration(tmpss.str(),SetSharedPtr(1,which_explaining_genotype));
						tmpss.str("");
						tmpss<<which_index<<"_"<<k<<"combinations";
						int f5=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,which_genotype));
						int f6=variables.getEnumeration("selections",SetSharedPtr(3,which_chromosome,which_explaining_genotype,which_genotype));
						res<<(-1*f4)<<" "<<f5<<" "<<(-1*f6)<<EOL;

					}
				}
			}
		}
	}

	for(int which_genotype=0;which_genotype<number_of_genotypes;which_genotype++){
		for(int which_explaining_genotype=0;which_explaining_genotype<number_of_explaining_haplotypes;which_explaining_genotype++){
			if(0 == which_explaining_genotype){
				for(int which_chromosome=0;which_chromosome<ploidy;which_chromosome++){
					int f1=variables.getEnumeration("v",SetSharedPtr(3,which_chromosome,which_explaining_genotype,which_genotype));
					int f2=variables.getEnumeration("selections",SetSharedPtr(3,which_chromosome,which_explaining_genotype,which_genotype));
					res<<(-1*f1)<<" "<<f2<<EOL;
					res<<f1<<" "<<(-1*f2)<<EOL;
				}
			}
			if(which_explaining_genotype<number_of_explaining_haplotypes-1){
				for(int which_chromosome=0;which_chromosome<ploidy;which_chromosome++){
					int f1=variables.getEnumeration("v",SetSharedPtr(3,which_chromosome,which_explaining_genotype+1,which_genotype));
					int f2=variables.getEnumeration("selections",SetSharedPtr(3,which_chromosome,which_explaining_genotype+1,which_genotype));
					int f3=variables.getEnumeration("v",SetSharedPtr(3,which_chromosome,which_explaining_genotype,which_genotype));
					res<<(-1*f1)<<" "<<f2<<" "<<f3<<EOL;
					int f4=variables.getEnumeration("selections",SetSharedPtr(3,which_chromosome,which_explaining_genotype+1,which_genotype));
					int f5=variables.getEnumeration("v",SetSharedPtr(3,which_chromosome,which_explaining_genotype+1,which_genotype));
					res<<(-1*f4)<<" "<<f5<<EOL;
					int f6=variables.getEnumeration("v",SetSharedPtr(3,which_chromosome,which_explaining_genotype,which_genotype));
					int f7=variables.getEnumeration("v",SetSharedPtr(3,which_chromosome,which_explaining_genotype+1,which_genotype));
					res<<(-1*f6)<<" "<<f7<<EOL;
				}
			}
			if(which_explaining_genotype>0 && (which_explaining_genotype<number_of_explaining_haplotypes-1)){
				for(int which_chromosome=0;which_chromosome<ploidy;which_chromosome++){
					int f1=variables.getEnumeration("v",SetSharedPtr(3,which_chromosome,which_explaining_genotype,which_genotype));
					int f2=variables.getEnumeration("selections",SetSharedPtr(3,which_chromosome,which_explaining_genotype+1,which_genotype));
					res<<(-1*f1)<<" "<<(-1*f2)<<EOL;
				}
			}
			if(which_explaining_genotype==number_of_explaining_haplotypes-1){
				for(int which_chromosome=0;which_chromosome<ploidy;which_chromosome++){
					int f1=variables.getEnumeration("v",SetSharedPtr(3,which_chromosome,which_explaining_genotype,which_genotype));
					res<<f1<<EOL;
				}
			}
		}
	}

	for(int which_explaining_haplotype=number_of_known_haplotypes;which_explaining_haplotype<number_of_explaining_haplotypes-1;which_explaining_haplotype++){
		int _index_e=0;
		for(int which_index=0;which_index<length_of_genotypes;which_index++){
			int width=CEIL(log2(this->occurence[0][which_index].size()));
			for(int k=0;k<width;k++){
				tmpss.str("");
				tmpss<<which_index<<"_"<<k<<"e";
				int f1=variables.getEnumeration(tmpss.str(),SetSharedPtr(1,which_explaining_haplotype));
				tmpss.str("");
				tmpss<<which_index<<"_"<<k<<"haplotypes";
				int f2=variables.getEnumeration(tmpss.str(),SetSharedPtr(1,which_explaining_haplotype));
				tmpss.str("");
				tmpss<<which_index<<"_"<<k<<"haplotypes";
				int f3=variables.getEnumeration(tmpss.str(),SetSharedPtr(1,which_explaining_haplotype+1));
				res<<f1<<" "<<(-1*f2)<<" "<<f3<<EOL;
				if(0 == which_index && 0 == k){
					tmpss.str("");
					tmpss<<which_index<<"_"<<k<<"e";
					int f1=variables.getEnumeration(tmpss.str(),SetSharedPtr(1,which_explaining_haplotype));
					tmpss.str("");
					tmpss<<which_index<<"_"<<k<<"haplotypes";
					int f2=variables.getEnumeration(tmpss.str(),SetSharedPtr(1,which_explaining_haplotype));
					res<<(-1*f1)<<" "<<(-1*f2)<<EOL;
					int f3=variables.getEnumeration(tmpss.str(),SetSharedPtr(1,which_explaining_haplotype+1));
					res<<(-1*f1)<<" "<<f3<<EOL;
					res<<f1<<" "<<f2<<" "<<(-1*f3)<<EOL;
				}
				if(which_index>0 || (0 == which_index && k>0)){//ln461
					tmpss.str("");
					tmpss<<which_index<<"_"<<k<<"e";
					int f1=variables.getEnumeration(tmpss.str(),SetSharedPtr(1,which_explaining_haplotype));
					tmpss.str("");
					tmpss<<which_index<<"_"<<k<<"haplotypes";
					int f2=variables.getEnumeration(tmpss.str(),SetSharedPtr(1,which_explaining_haplotype));
					res<<(-1*f1)<<" "<<_index_e<<" "<<(-1*f2)<<EOL;
					int f3=variables.getEnumeration(tmpss.str(),SetSharedPtr(1,which_explaining_haplotype+1));
					res<<(-1*f1)<<" "<<_index_e<<" "<<f3<<EOL;
					res<<f1<<" "<<f2<<" "<<(-1*f3)<<EOL;
					res<<f1<<" "<<(-1*_index_e)<<EOL;
				}
				tmpss.str("");
				tmpss<<which_index<<"_"<<k<<"e";
				int f=variables.getEnumeration(tmpss.str(),SetSharedPtr(1,which_explaining_haplotype));
				if((which_index == length_of_genotypes-1) && (width-1 == k)){
					res<<f<<EOL;
				}
					_index_e=f;


			}

		}
	}

	//ln497
	for(int which_genotype=0; which_genotype<number_of_genotypes;which_genotype++){
		for(int which_chromosome=0;which_chromosome<ploidy-1;which_chromosome++){
			int _index_s=0;
			for(int which_explaining_haplotype=0;which_explaining_haplotype<number_of_explaining_haplotypes;which_explaining_haplotype++){
				int f1=variables.getEnumeration("s",SetSharedPtr(3,which_chromosome,which_explaining_haplotype,which_genotype));
				int f2=variables.getEnumeration("selections",SetSharedPtr(3,which_chromosome,which_explaining_haplotype,which_genotype));
				int f3=variables.getEnumeration("selections",SetSharedPtr(3,which_chromosome+1,which_explaining_haplotype,which_genotype));
				res<<f1<<" "<<(-1*f2)<<" "<<f3<<EOL;
				if(0 == which_explaining_haplotype){
					int f1=variables.getEnumeration("s",SetSharedPtr(3,which_chromosome,which_explaining_haplotype,which_genotype));
					int f2=variables.getEnumeration("selections",SetSharedPtr(3,which_chromosome,which_explaining_haplotype,which_genotype));
					res<<(-1*f1)<<" "<<(-1*f2)<<EOL;
					int f3=variables.getEnumeration("selections",SetSharedPtr(3,which_chromosome+1,which_explaining_haplotype,which_genotype));
					res<<(-1*f1)<<" "<<f3<<EOL;
					res<<f1<<" "<<f2<<" "<<(-1*f3)<<EOL;
				};
				if(which_explaining_haplotype>0){
					int f1=variables.getEnumeration("s",SetSharedPtr(3,which_chromosome,which_explaining_haplotype,which_genotype));
					int f2=variables.getEnumeration("selections",SetSharedPtr(3,which_chromosome,which_explaining_haplotype,which_genotype));
					res<<(-1*f1)<<" "<<_index_s<<" "<<(-1*f2)<<EOL;
					int f3=variables.getEnumeration("selections",SetSharedPtr(3,which_chromosome+1,which_explaining_haplotype,which_genotype));
					res<<(-1*f1)<<" "<<_index_s<<" "<<f3<<EOL;
					res<<f1<<" "<<f2<<" "<<(-1*f3)<<EOL;
					res<<f1<<" "<<(-1*_index_s)<<EOL;
				}
				_index_s=variables.getEnumeration("s",SetSharedPtr(3,which_chromosome,which_explaining_haplotype,which_genotype));
			}
			//ln559
			for(int which_explaining_haplotype=0;which_explaining_haplotype<number_of_explaining_haplotypes;which_explaining_haplotype++){
				int f1=variables.getEnumeration("selections",SetSharedPtr(3,which_chromosome,which_explaining_haplotype,which_genotype));
				int f2=variables.getEnumeration("selections",SetSharedPtr(3,which_chromosome+1,which_explaining_haplotype,which_genotype));
				res<<_index_s<<" "<<(-1*f1)<<" "<<f2<<EOL;
				res<<_index_s<<" "<<f1<<" "<<(-1*f2)<<EOL;
			}
		}
	}
//	this->VarNum=variables.getVarnum();
//	std::cout<<"\n"<<res.str();
//	std::cout<<"\n"<<this->ClauseNum<<","<<this->VarNum<<"\n";

};

void Haplotype::getBiallelicCoding(int ploidy,int which_genotype, int which_index, int which_allele, IndexingVariables& variables){
	int allele=this->occurence[which_genotype][which_index][which_allele];
	std::stringstream tmpss;
	int q=MIN(allele,ploidy-allele);
	int q_=CEIL(log2(q+1));
	q_=q_<1?1:q_;
	for(int t=0;t<q_;t++){
		tmpss.str("");
		tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"sum";
		int a=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,0,t));
		res<<(-1*a)<<EOL;
	}
	for(int l=0;l<ploidy;l++){
		tmpss.str("");
		tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
		int c=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,l,q_));
		res<<(-1*c)<<EOL;
	}
	std::string s=ToBinaryString(q);
	boost::shared_ptr<bool[]> b(new bool[q_]);
	for(int i=0;i<q_;i++){
		b[i]=false;
	}
	int index_=0;
	for(int index=s.length()-1;index>=0;index--){
		if(s[index]=='1')
			b[index_]=true;
		index_++;
	}
	for(int t=0;t<q_;t++)
	{
		tmpss.str("");
		tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"sum";
		int a=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,ploidy,t));
		if(b[t]){
			res<<a<<EOL;}
		else{
			res<<(-1)*a<<EOL;}
	}
	for(int l=0;l<ploidy;l++){
		for(int t=0;t<q_;t++){
			tmpss.str("");
			tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"sum";
			int a_=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,l+1,t));
			int a =variables.getEnumeration(tmpss.str(),SetSharedPtr(2,l,t));
			tmpss.str("");
			tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
			int c=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,l,t));
			int c_=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,l,t+1));

			res<<a_<<" "<<a<<" "<<(-1*c)<<EOL;
			res<<a_<<" "<<(-1*a)<<" "<<c<<EOL;
			res<<(-1*a_)<<" "<<(-1*a)<<" "<<(-1*c)<<EOL;
			res<<(-1*a_)<<" "<<a<<" "<<c<<EOL;
			res<<c_<<" "<<(-1*a)<<" "<<(-1*c)<<EOL;
			res<<(-1*c_)<<" "<<a<<EOL;
			res<<(-1*c_)<<" "<<c<<EOL;
		}
	}
}

void Haplotype::getGeneralCoding(int ploidy, int which_genotype, int which_index, IndexingVariables& variables){
	std::stringstream tmpss;
	int number_of_different_alleles=this->occurence[which_genotype][which_index].size();
	int width=CEIL(log2(this->occurence[which_genotype][which_index].size()));
	std::vector< boost::shared_ptr< int[]> > mbool(number_of_different_alleles);

	for(int which_allele=0;which_allele<number_of_different_alleles;which_allele++){
		mbool[which_allele]=toBooleanInt(width,which_allele);
	}
	for(int which_allele=0;which_allele<number_of_different_alleles-1;which_allele++){
		int min=MIN(this->occurence[which_genotype][which_index][which_allele],ploidy-this->occurence[which_genotype][which_index][which_allele]);
		if(this->occurence[which_genotype][which_index][which_allele] == min){
			for(int which_chromosome=0;which_chromosome<ploidy;which_chromosome++){

				for(int k=0;k<width;k++){
					int f1=(int)pow(-1,mbool[which_allele][k]);
					tmpss.str("");;
					tmpss<<which_index<<"_"<<k<<"combinations";
					int f2=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,which_genotype));
					res<<(f1*f2)<<" ";
				}
				tmpss.str("");;
				tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
					res<<variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0))<<EOL;
				for(int k=0;k<width;k++){
						int f1=pow(-1,1-mbool[which_allele][k]);
					tmpss.str("");;
					tmpss<<which_index<<"_"<<k<<"combinations";
					int f2=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,which_genotype));
					res<<(f1*f2)<<" ";
					tmpss.str("");;
					tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
					int f3=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
					res<<(-1*f3)<<EOL;
				}
			}
		}else{
			for(int which_chromosome=0;which_chromosome<ploidy;which_chromosome++){
				for(int k=0;k<width;k++){
					tmpss.str("");;
					tmpss<<which_index<<"_"<<k<<"combinations";
					int f1=pow(-1,mbool[which_allele][k]);
					int f2=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,which_genotype));
					res<<f1*f2<<" ";
				}
				tmpss.str("");;
				tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
				int f1=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
				res<<(-1)*f1<<EOL;
				for(int k=0;k<width;k++){
					int f1=pow(-1,1-mbool[which_allele][k]);
					tmpss.str("");;
					tmpss<<which_index<<"_"<<k<<"combinations";
					int f2=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,which_genotype));
					res<<(f1*f2)<<" ";
					tmpss.str("");;
					tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
					int f3=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
					res<<f3<<EOL;
				}
			}
		}
	}
		for(int which_chromosome=0;which_chromosome<ploidy;which_chromosome++){
			for(int k=0;k<width;k++){
				int f1=pow(-1,mbool[number_of_different_alleles-1][k]);
				tmpss.str("");;
				tmpss<<which_index<<"_"<<k<<"combinations";
				int f2=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,which_genotype));
				res<<(f1*f2)<<" ";
			}
			tmpss.str("");;
			tmpss<<which_genotype<<"_"<<which_index<<"_"<<(number_of_different_alleles-1)<<"parity";
			int f=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
			res<<f<<EOL;
			for(int k=0;k<width;k++){
				tmpss.str("");;
				int f1=pow(-1,1-mbool[number_of_different_alleles-1][k]);
				tmpss<<which_index<<"_"<<k<<"combinations";
				int f2=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,which_genotype));
				tmpss.str("");;
				tmpss<<which_genotype<<"_"<<which_index<<"_"<<(number_of_different_alleles-1)<<"parity";
				int f3=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
				res<<(f1*f2)<<" "<<(-1)*f3<<EOL;
			}
		}

		for(int which_chromosome=0;which_chromosome<ploidy;which_chromosome++){
			for(int which_allele=0;which_allele<number_of_different_alleles-1;which_allele++){
				int min=MIN(this->occurence[which_genotype][which_index][which_allele],ploidy-this->occurence[which_genotype][which_index][which_allele]);
				if(min==this->occurence[which_genotype][which_index][which_allele]){
					tmpss.str("");;
					tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
					int f=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
					res<<f<<" ";
				}else{
					tmpss.str("");;
					tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
					int f=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
					res<<(-1*f)<<" ";
				}
			};
				tmpss.str("");;
				tmpss<<which_genotype<<"_"<<which_index<<"_"<<(number_of_different_alleles-1)<<"parity";
				int f=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
				res<<f<<EOL;
				for(int which_allele=0;which_allele<number_of_different_alleles-1;which_allele++){
					int min=MIN(this->occurence[which_genotype][which_index][which_allele],ploidy-this->occurence[which_genotype][which_index][which_allele]);
					if(min==this->occurence[which_genotype][which_index][which_allele]){
						tmpss.str("");;
						tmpss<<which_genotype<<"_"<<which_index<<"_"<<(number_of_different_alleles-1)<<"parity";
						int f1=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
						tmpss.str("");;
						tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
						int f2=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
						res<<(-1*f1)<<" "<<(-1*f2)<<EOL;
					}else{
						tmpss.str("");;
						tmpss<<which_genotype<<"_"<<which_index<<"_"<<(number_of_different_alleles-1)<<"parity";
						int f1=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
						tmpss.str("");;
						tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
						int f2=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
						res<<(-1*f1)<<" "<<f2<<EOL;

				}

			}
		}

		for(int which_allele=0;which_allele<number_of_different_alleles-1;which_allele++){
			getBiallelicCoding(ploidy,which_genotype,which_index,which_allele,variables);
		}

}

void Haplotype::getGeneralCodingMissing(int ploidy, int which_genotype, int which_index, IndexingVariables& variables){
	std::stringstream tmpss;
	int number_of_different_alleles=this->occurence[which_genotype][which_index].size();
	int width=CEIL(log2(this->occurence[which_genotype][which_index].size()));
	std::vector< boost::shared_ptr< int[]> > mbool(number_of_different_alleles);
	for(int which_allele=0;which_allele<number_of_different_alleles;which_allele++){
		mbool[which_allele]=toBooleanInt(width,which_allele);
	}

	for(int which_allele=0;which_allele<number_of_different_alleles;which_allele++){
		for(int which_chromosome=0;which_chromosome<ploidy;which_chromosome++){
			for(int k=0;k<width;k++){
				int f1=pow(-1,mbool[which_allele][k]);
				tmpss.str("");
				tmpss<<which_index<<"_"<<k<<"combinations";
				int f2=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,which_genotype));
				res<<f1*f2<<" ";
			}
			tmpss.str("");
			tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
			int f1=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
			res<<f1<<EOL;
			for(int k=0;k<width;k++){
				int f1=pow(-1,1-mbool[which_allele][k]);
				tmpss.str("");
				tmpss<<which_index<<"_"<<k<<"combinations";
				int f2=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,which_genotype));
				tmpss.str("");
				tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
				int f3=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
				res<<f1*f2<<" "<<(-1*f3)<<EOL;
			}
		}
	}
//ln206
	for(int which_allele=0;which_allele<number_of_different_alleles;which_allele++){
		if(this->occurence[which_genotype][which_index][which_allele]>0){
			for(int which_chromosome=0;which_chromosome<ploidy;which_chromosome++){
				tmpss.str("");
				tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
				int f1=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
				res<<f1<<" ";
			}
			res<<EOL2;
		}else{
			for(int which_chromosome=0;which_chromosome<ploidy;which_chromosome++){
				tmpss.str("");
				tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
				int f1=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
				res<<(-1*f1)<<EOL;
			}
		}
	}

	//ln233
	for(int which_chromosome=0;which_chromosome<ploidy;which_chromosome++){
		for(int which_allele=0;which_allele<number_of_different_alleles;which_allele++){
			tmpss.str("");
			tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
			int f1=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
			res<<f1<<" ";
		}
		res<<EOL2;
	}
}


void Haplotype::getGeneralCodingTotalyMissing(int ploidy, int which_genotype, int which_index, IndexingVariables& variables){
	std::stringstream tmpss;
	int number_of_different_alleles=this->occurence[which_genotype][which_index].size();
	int width=CEIL(log2(this->occurence[which_genotype][which_index].size()));
	std::vector< boost::shared_ptr< int[]> > mbool(number_of_different_alleles);
	for(int which_allele=0;which_allele<number_of_different_alleles;which_allele++){
		mbool[which_allele]=toBooleanInt(width,which_allele);
	}

	//262
	for(int which_allele=0;which_allele<number_of_different_alleles;which_allele++){
		for(int which_chromosome=0;which_chromosome<ploidy;which_chromosome++){
			for(int k=0;k<width;k++){
				int f1=pow(-1,mbool[which_allele][k]);
				tmpss.str("");
				tmpss<<which_index<<"_"<<k<<"combinations";
				int f2=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,which_genotype));
				res<<f1*f2<<" ";
			}
			tmpss.str("");
			tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
			int f1=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
			res<<f1<<EOL;
			for(int k=0;k<width;k++){
				int f1=pow(-1,1-mbool[which_allele][k]);
				tmpss.str("");
				tmpss<<which_index<<"_"<<k<<"combinations";
				int f2=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,which_genotype));
				tmpss.str("");
				tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
				int f3=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
				res<<f1*f2<<" "<<(-1*f3)<<EOL;
			}
		}
	}

	//ln290
	for(int which_chromosome=0;which_chromosome<ploidy;which_chromosome++){
		for(int which_allele=0;which_allele<number_of_different_alleles;which_allele++){
			tmpss.str("");
			tmpss<<which_genotype<<"_"<<which_index<<"_"<<which_allele<<"parity";
			int f1=variables.getEnumeration(tmpss.str(),SetSharedPtr(2,which_chromosome,0));
			res<<f1<<" ";
		}
		res<<EOL2;
	};

}


void Haplotype::createVariables(int number_of_explaining_haplotypes,IndexingVariables& variables){
	int number_of_genotypes=this->occurence.shape()[0];
	int length_of_genotypes=this->occurence.shape()[1];
	int ploidy=this->data->getNumOfChrSet();
	boost::shared_ptr<int[]> sp;
	sp=SetSharedPtr(3,ploidy,number_of_explaining_haplotypes,number_of_genotypes);
	variables.add("selections",sp);
	variables.add("s",sp);
	variables.add("v",sp);
	std::stringstream tmpss;
	for(int which_index=0;which_index<length_of_genotypes;which_index++){
		int width=CEIL(log2(this->occurence[0][which_index].size()));
		for(int k=0;k<width;k++){
			tmpss.str("");
			tmpss<<which_index<<"_"<<k<<"haplotypes";
			variables.add(tmpss.str(),SetSharedPtr(1,number_of_explaining_haplotypes));
			tmpss.str("");
			tmpss<<which_index<<"_"<<k<<"combinations";
			variables.add(tmpss.str(),SetSharedPtr(2,ploidy,number_of_genotypes));
			tmpss.str("");
			tmpss<<which_index<<"_"<<k<<"e";
			variables.add(tmpss.str(),SetSharedPtr(1,number_of_explaining_haplotypes));
		};
	};

	int counter=0;
	for(int which_index=0;which_index<length_of_genotypes;which_index++){
		for(int which_allele=0;which_allele<this->occurence[0][which_index].size()-1;which_allele++){
			tmpss.str("");;
			tmpss<<which_index<<"_"<<which_allele<<"h";
			variables.add(tmpss.str(),SetSharedPtr(1,number_of_explaining_haplotypes));
			counter++;
		}
	}

	int c=(int)(0.5*counter*(counter-1));
	sp=SetSharedPtr(2,number_of_explaining_haplotypes,c);
	variables.add("a",sp);
	variables.add("b",sp);
	variables.add("c",sp);

	sp=SetSharedPtr(1,c);
	variables.add("A",sp);
	variables.add("B",sp);
	variables.add("C",sp);
	variables.add("H",sp);

	int Q=-1;
	int Q_=CEIL(log2(Q+1));
	if(Q_<1)
		Q_=1;
	variables.add("bla",SetSharedPtr(2,c+1,Q_));
	variables.add("blu",SetSharedPtr(2,c,Q_+2));

	for(int which_genotype=0;which_genotype<number_of_genotypes;which_genotype++){
		for(int which_index=0;which_index<length_of_genotypes;which_index++){
			for(int k=0;k<this->occurence[which_genotype][which_index].size();k++){
				if(0 == missing[which_genotype][which_index]){
					int allele=this->occurence[which_genotype][which_index][k];
					int q = MIN(allele,ploidy-allele);
					int q_=CEIL(log2(q+1));
					if(q_<1)
						q_=1;
					tmpss.str("");;
					tmpss<<which_genotype<<"_"<<which_index<<"_"<<k<<"sum";
					variables.add(tmpss.str(),SetSharedPtr(2,ploidy+1,q_));
					tmpss.str("");
					tmpss<<which_genotype<<"_"<<which_index<<"_"<<k<<"parity";
					variables.add(tmpss.str(),SetSharedPtr(2,ploidy,q_+2));
				}else
				{
					tmpss.str("");
					tmpss<<which_genotype<<"_"<<which_index<<"_"<<k<<"parity";
					variables.add(tmpss.str(),SetSharedPtr(2,ploidy,1));
				}
			}
		}
	}
}

} /* namespace SHEsis */
