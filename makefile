
#LIB= /home/ada/git/SHEsis/boost/boost_1_55_0/stage/lib
#INC= /home/ada/git/SHEsis/boost/boost_1_55_0/
#CXXFLAGS = -O0 -g  -std=c++11 -fmessage-length=0   -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS  -Wno-parentheses -L$(LIB) -I$(INC)
LIB= -L/results/software/boost_1_55_0/stage/lib -L/results/software/mlpack-1.0.11/built/lib -L/results/software/armadillo-4.500.1
INC= -I/results/software/boost_1_55_0/ -I/results/software/mlpack-1.0.11/built/include -I/results/software/armadillo-4.500.1/include -I/results/software/libxml2-2.9.2/include
CXXFLAGS =  -fmessage-length=0   -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS  -Wno-parentheses $(LIB) $(INC) #-pg 
 ifdef RELEASE
 CXXFLAGS += -O3 
 else
 CXXFLAGS += -O0 -g
 endif

TARGET = SHEsis SHEsisData_test AssociationTest_test HWETest_test  LDTest_test QTL_test DataGenerator GeneInteractionQTL_test GeneInteractionBinary_test logistic linear MarkerRegression_test HaplotypeEM_test HaplotypeLD_test

all: QTL_test HaplotypeDiploid_test SHEsis SHEsisData_test AssociationTest_test Haplotype_test HWETest_test LDTest_test HaplotypeEM_test DataGenerator logistic linear MarkerRegression_test  HaplotypeLD_test


SHEsis:	 main.o  HaplotypeLD.o SHEsisData.o fisher.o utility.o AssociationTest.o HWETest.o LDTest.o QTL.o HaplotypeBase.o Haplotype.o HaplotypeEM.o  IndexingVariables.o  ArrayStorage.o System.o Solver.o Options.o  BMP.o font.o minifont.o CreatHtmlTable.o linear.o regression.o logistic.o  MarkerRegression.o GeneInteractionQTL.o GeneInteractionBinary.o GeneInteraction.o
	$(CXX) $(LIB) -lmlpack   -o SHEsis SHEsisData.o HaplotypeLD.o main.o fisher.o utility.o AssociationTest.o HWETest.o LDTest.o HaplotypeBase.o HaplotypeEM.o Haplotype.o CreatHtmlTable.o  IndexingVariables.o ArrayStorage.o System.o Solver.o Options.o  BMP.o font.o minifont.o QTL.o linear.o regression.o logistic.o  MarkerRegression.o  GeneInteractionQTL.o GeneInteractionBinary.o GeneInteraction.o

DataGenerator: DataGenerator.o
	$(CXX) $(LIB) -o DataGenerator DataGenerator.o -l:libboost_program_options.a 

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

HaplotypeEM_test: HaplotypeEM.o SHEsisData.o HaplotypeBase.o HaplotypeEM_test.o fisher.o utility.o CreatHtmlTable.o
	$(CXX) -o HaplotypeEM_test HaplotypeEM.o SHEsisData.o HaplotypeBase.o HaplotypeEM_test.o fisher.o utility.o CreatHtmlTable.o

HaplotypeLD_test: HaplotypeLD_test.o BMP.o font.o minifont.o  Haplotype.o IndexingVariables.o ArrayStorage.o System.o Solver.o Options.o   LDTest.o HaplotypeLD.o HaplotypeEM.o SHEsisData.o HaplotypeBase.o  fisher.o utility.o CreatHtmlTable.o
	$(CXX) -o HaplotypeLD_test BMP.o HaplotypeLD_test.o font.o minifont.o  Haplotype.o IndexingVariables.o ArrayStorage.o System.o Solver.o Options.o   LDTest.o  HaplotypeLD.o HaplotypeEM.o SHEsisData.o HaplotypeBase.o fisher.o utility.o CreatHtmlTable.o

LDTest_test: SHEsisData.o LDTest_test.o HaplotypeEM.o LDTest.o HaplotypeBase.o BMP.o font.o minifont.o Haplotype.o  IndexingVariables.o ArrayStorage.o System.o Solver.o Options.o  fisher.o utility.o  CreatHtmlTable.o
	$(CXX) -o LDTest_test SHEsisData.o HaplotypeBase.o HaplotypeEM.o LDTest_test.o LDTest.o BMP.o font.o minifont.o Haplotype.o IndexingVariables.o ArrayStorage.o System.o Solver.o Options.o   fisher.o utility.o  CreatHtmlTable.o

QTL_test: SHEsisData.o  utility.o QTL.o QTL_test.o CreatHtmlTable.o
	$(CXX) -o QTL_test SHEsisData.o utility.o QTL.o QTL_test.o CreatHtmlTable.o

logistic: logistic.o logistic_test.o regression.o 
	$(CXX)  $(LIB) -lmlpack  -o logistic_test logistic.o logistic_test.o regression.o

linear: linear.o linear_test.o regression.o 
	$(CXX)  $(LIB) -lmlpack  -o linear_test linear.o linear_test.o regression.o

MarkerRegression_test: linear.o regression.o logistic.o SHEsisData.o MarkerRegression_test.o MarkerRegression.o utility.o CreatHtmlTable.o
	$(CXX)  $(LIB) -lmlpack  -o MarkerRegression_test linear.o regression.o logistic.o SHEsisData.o MarkerRegression_test.o MarkerRegression.o utility.o CreatHtmlTable.o

GeneInteractionBinary_test:  GeneInteractionBinary.o SHEsisData.o GeneInteraction.o GeneInteractionBinary_test.o utility.o CreatHtmlTable.o
	$(CXX)   $(LIB)  -o GeneInteractionBinary_test  GeneInteractionBinary.o SHEsisData.o GeneInteraction.o GeneInteractionBinary_test.o utility.o CreatHtmlTable.o

GeneInteractionQTL_test: SHEsisData.o GeneInteraction.o utility.o linear.o regression.o GeneInteractionQTL.o GeneInteractionQTL_test.o CreatHtmlTable.o
	$(CXX)  $(LIB) -lmlpack  -o GeneInteractionQTL_test SHEsisData.o GeneInteraction.o utility.o linear.o regression.o GeneInteractionQTL.o GeneInteractionQTL_test.o CreatHtmlTable.o

clean:
	rm -f *.o $(TARGET)
