CXXFLAGS =	-O2 -g -Wall -fmessage-length=0  -L/home/ada/git/SHEsis/boost/boost_1_55_0/stage/lib -I/home/ada/git/SHEsis/boost/boost_1_55_0 

LIBS =


SHEsis:	SHEsisData.o main.o
	$(CXX) -o SHEsis SHEsisData.o main.o $(LIBS)

 
all:	SHEsis SHEsisData_test

SHEsisData_test: SHEsisData.o SHEsisData_test.o
	$(CXX) -o  SHEsisData_test SHEsisData.o SHEsisData_test.o $(LIBS)

clean:
	rm -f $(OBJS) $(TARGET)
