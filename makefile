CXXFLAGS =	 -O0  -g -Wall -std=c++11 -fmessage-length=0  -I/home/ada/git/SHEsis/boost/boost_1_55_0 -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS  -Wno-parentheses -L/home/ada/git/SHEsis/boost/boost_1_55_0/stage/lib 
#CXXFLAGS =	-O0 -g -Wall -fmessage-length=0   -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS  -Wno-parentheses -L/results/software/boost_1_55_0/stage/lib -I/results/software/boost_1_55_0/ 
 
LIBS = 
TARGET = SHEsis SHEsisData_test AssociationTest_test HWETest_test HaplotypeDiploid_test
OBJS = SHEsisData.o main.o SHEsisData_test.o utility.o fisher.o AssociationTest.o AssociationTest_test.o \
HWETest.o  HWETest_test.o Haplotype.o 


all:	HaplotypeDiploid_test SHEsis SHEsisData_test AssociationTest_test Haplotype_test HWETest_test

SHEsis:	 main.o  SHEsisData.o fisher.o utility.o AssociationTest.o HWETest.o LDTest.o Haplotype.o HaplotypeDiploid.o IndexingVariables.o ArrayStorage.o System.o Solver.o Options.o
	$(CXX) -L/home/ada/git/SHEsis/boost/boost_1_55_0/stage/lib  -o SHEsis SHEsisData.o main.o fisher.o utility.o AssociationTest.o HWETest.o LDTest.o Haplotype.o HaplotypeDiploid.o \
	 IndexingVariables.o ArrayStorage.o System.o Solver.o Options.o -l:libboost_program_options.a

SHEsisData_test: SHEsisData.o SHEsisData_test.o utility.o
	$(CXX) -o  SHEsisData_test SHEsisData.o SHEsisData_test.o $(LIBS)

AssociationTest_test: SHEsisData.o utility.o fisher.o AssociationTest.o AssociationTest_test.o 
	$(CXX) -o AssociationTest_test SHEsisData.o utility.o fisher.o AssociationTest.o AssociationTest_test.o $(LIBS)

HWETest_test: SHEsisData.o utility.o fisher.o HWETest.o HWETest_test.o 
	$(CXX) -o HWETest_test SHEsisData.o utility.o fisher.o HWETest.o HWETest_test.o $(LIBS)

Haplotype_test: SHEsisData.o    Haplotype_test.o  fisher.o utility.o Haplotype.o IndexingVariables.o ArrayStorage.o System.o Solver.o Options.o
	$(CXX) -o Haplotype_test SHEsisData.o    Haplotype_test.o  fisher.o utility.o Haplotype.o IndexingVariables.o ArrayStorage.o System.o Solver.o Options.o $(LIBS)


HaplotypeDiploid_test: SHEsisData.o HaplotypeDiploid_test.o fisher.o utility.o HaplotypeDiploid.o 
	$(CXX) -o HaplotypeDiploid_test SHEsisData.o HaplotypeDiploid_test.o fisher.o utility.o HaplotypeDiploid.o  $(LIBS)


clean:
	rm -f *.o $(TARGET)
