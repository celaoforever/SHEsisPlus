
LIB= /path/to/boost_1_55_0/stage/lib
INC= /path/to/software/boost_1_55_0/ 
CXXFLAGS = -O3  -fmessage-length=0   -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS  -Wno-parentheses -L$(LIB) -I$(INC)

TARGET = SHEsisPlus 

all: $(TARGET)

SHEsisPlus:	 main.o  SHEsisData.o fisher.o utility.o AssociationTest.o HWETest.o LDTest.o QTL.o HaplotypeBase.o Haplotype.o HaplotypeEM.o  IndexingVariables.o ArrayStorage.o System.o Solver.o Options.o  BMP.o font.o minifont.o CreatHtmlTable.o 
	$(CXX) -L$(LIB)  -o SHEsisPlus SHEsisData.o main.o fisher.o utility.o AssociationTest.o HWETest.o LDTest.o HaplotypeBase.o HaplotypeEM.o Haplotype.o CreatHtmlTable.o\
	 IndexingVariables.o ArrayStorage.o System.o Solver.o Options.o  BMP.o font.o minifont.o QTL.o -l:libboost_program_options.a 

clean:
	rm -f *.o $(TARGET)
