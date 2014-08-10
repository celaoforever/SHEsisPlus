CXXFLAGS =	-O2 -g -Wall -fmessage-length=0  -L/results/software/boost_1_55_0/stage/lib -I/results/software/boost_1_55_0/ 

LIBS =
TARGET = SHEsis SHEsisData_test AssociationTest_test
OBJS = SHEsisData.o main.o SHEsisData_test.o utility.o fisher.o AssociationTest.o AssociationTest_test.o


all:	SHEsis SHEsisData_test AssociationTest_test

SHEsis:	 main.o  SHEsisData.o fisher.o utility.o AssociationTest.o
	$(CXX) -o SHEsis SHEsisData.o main.o $(LIBS)

SHEsisData_test: SHEsisData.o SHEsisData_test.o utility.o
	$(CXX) -o  SHEsisData_test SHEsisData.o SHEsisData_test.o $(LIBS)

AssociationTest_test: SHEsisData.o utility.o fisher.o AssociationTest.o AssociationTest_test.o 
	$(CXX) -o AssociationTest_test SHEsisData.o utility.o fisher.o AssociationTest.o AssociationTest_test.o $(LIBS)

clean:
	rm -f $(OBJS) $(TARGET)
