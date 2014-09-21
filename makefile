
LIB= /results/software/boost_1_55_0/stage/lib
INC= /results/software/boost_1_55_0/ 

#CXXFLAGS =	 -O0  -g -Wall -std=c++11 -fmessage-length=0  -I/home/ada/git/SHEsis/boost/boost_1_55_0 -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS  -Wno-parentheses -L/home/ada/git/SHEsis/boost/boost_1_55_0/stage/lib 
CXXFLAGS =	-O0 -g -Wall -fmessage-length=0   -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS  -Wno-parentheses -L$(LIB) -I$(INC)
 

TARGET = SHEsis SHEsisData_test AssociationTest_test HWETest_test HaplotypeDiploid_test LDTest_test
OBJS = SHEsisData.o main.o SHEsisData_test.o utility.o fisher.o AssociationTest.o AssociationTest_test.o \
HWETest.o  HWETest_test.o Haplotype.o  CreatHtmlTable.o


all:	HaplotypeDiploid_test SHEsis SHEsisData_test AssociationTest_test Haplotype_test HWETest_test LDTest_test

SHEsis:	 main.o  SHEsisData.o fisher.o utility.o AssociationTest.o HWETest.o LDTest.o HaplotypeBase.o Haplotype.o HaplotypeDiploid.o IndexingVariables.o ArrayStorage.o System.o Solver.o Options.o  BMP.o font.o minifont.o CreatHtmlTable.o 
	$(CXX) -L$(LIB)  -o SHEsis SHEsisData.o main.o fisher.o utility.o AssociationTest.o HWETest.o LDTest.o HaplotypeBase.o Haplotype.o HaplotypeDiploid.o CreatHtmlTable.o\
	 IndexingVariables.o ArrayStorage.o System.o Solver.o Options.o  BMP.o font.o minifont.o  -l:libboost_program_options.a 

SHEsisData_test: SHEsisData.o SHEsisData_test.o utility.o
	$(CXX) -o  SHEsisData_test SHEsisData.o SHEsisData_test.o $(LIBS)

AssociationTest_test: SHEsisData.o utility.o fisher.o AssociationTest.o AssociationTest_test.o  CreatHtmlTable.o
	$(CXX) -o AssociationTest_test SHEsisData.o utility.o fisher.o AssociationTest.o AssociationTest_test.o CreatHtmlTable.o $(LIBS)

HWETest_test: SHEsisData.o utility.o fisher.o HWETest.o HWETest_test.o  CreatHtmlTable.o
	$(CXX) -o HWETest_test SHEsisData.o utility.o fisher.o HWETest.o HWETest_test.o CreatHtmlTable.o $(LIBS)

Haplotype_test: SHEsisData.o  HaplotypeBase.o  Haplotype_test.o  fisher.o utility.o Haplotype.o IndexingVariables.o ArrayStorage.o System.o Solver.o Options.o CreatHtmlTable.o
	$(CXX) -o Haplotype_test SHEsisData.o  HaplotypeBase.o  Haplotype_test.o  fisher.o utility.o Haplotype.o IndexingVariables.o ArrayStorage.o System.o Solver.o Options.o CreatHtmlTable.o $(LIBS)


HaplotypeDiploid_test: SHEsisData.o HaplotypeBase.o HaplotypeDiploid_test.o fisher.o utility.o HaplotypeDiploid.o CreatHtmlTable.o
	$(CXX) -o HaplotypeDiploid_test SHEsisData.o HaplotypeBase.o HaplotypeDiploid_test.o fisher.o utility.o HaplotypeDiploid.o  CreatHtmlTable.o $(LIBS)
	
LDTest_test: SHEsisData.o LDTest_test.o LDTest.o HaplotypeBase.o BMP.o font.o minifont.o Haplotype.o HaplotypeDiploid.o IndexingVariables.o ArrayStorage.o System.o Solver.o Options.o  fisher.o utility.o  CreatHtmlTable.o
	$(CXX) -o LDTest_test SHEsisData.o HaplotypeBase.o LDTest_test.o LDTest.o BMP.o font.o minifont.o Haplotype.o HaplotypeDiploid.o IndexingVariables.o ArrayStorage.o System.o Solver.o Options.o   fisher.o utility.o  CreatHtmlTable.o

clean:
	rm -f *.o $(TARGET)
