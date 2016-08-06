CXX=g++
CXXFLAGS=-Wall -g -std=c++11 

all: Numerov-MB Numerov-1B+NB Numerov-1B Numerov-no-dip

Numerov-MB: main_MB.o spline.o wavefunction.o
	$(CXX) $(CXXFLAGS) -o Numerov-MB main_MB.o spline.o wavefunction.o

Numerov-1B+NB: main_1B+NB.o spline.o wavefunction.o
	$(CXX) $(CXXFLAGS) -o Numerov-1B+NB main_1B+NB.o spline.o wavefunction.o

Numerov-no-dip: main_no-dip.o spline.o wavefunction.o
	$(CXX) $(CXXFLAGS) -o Numerov-no-dip main_no-dip.o spline.o wavefunction.o

Numerov-1B: main_1B.o spline.o wavefunction.o
	$(CXX) $(CXXFLAGS) -o Numerov-1B main_1B.o spline.o wavefunction.o

main_1B+NB.o: main_1B+NB.cpp spline.h
	$(CXX) $(CXXFLAGS) -c main_1B+NB.cpp

main_MB.o: main_MB.cpp spline.h
	$(CXX) $(CXXFLAGS) -c main_MB.cpp

main_1B.o: main_1B.cpp spline.h
	$(CXX) $(CXXFLAGS) -c main_1B.cpp

main_no-dip.o: main_no-dip.cpp spline.h
	$(CXX) $(CXXFLAGS) -c main_no-dip.cpp

spline.o: spline.cpp
	$(CXX) $(CXXFLAGS) -c spline.cpp

wavefunction.o: wavefunction.cpp
	$(CXX) $(CXXFLAGS) -c wavefunction.cpp

clean: 
	rm *.o Numerov-1B Numerov-1B+NB Numerov-MB
