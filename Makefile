CXX=g++

all: Numerov-analysis

Numerov-analysis: main.o spline.o
	$(CXX) -o Numerov-analysis main.o spline.o

main.o: main.cpp spline.h
	$(CXX) -c main.cpp

spline.o: spline.cpp
	$(CXX) -c spline.cpp
