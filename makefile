CXXFLAGS =	-O2 -g -Wall -fmessage-length=0  -L/results/software/boost_1_55_0/stage/lib -I/results/software/boost_1_55_0/ 
 
OBJS = GenotypeMatrix.o /SHEsisData.o main.o
LIBS =

TARGET =	SHEsis

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)


all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
