IDIR = include
EIGDIR = ../Eigen
CXXFLAGS = -std=c++11

all : build/bin/openDGM

build/bin/openDGM : build/obj/functions1D.o build/obj/edge.o build/obj/element.o
	mkdir -p build/bin
	g++ $(CXXFLAGS) -I$(IDIR) -I$(EIGDIR) build/obj/functions1D.o build/obj/edge.o build/obj/element.o main.cpp -o $@

build/obj/edge.o : source/edge.cpp
	mkdir -p build/obj
	g++ -c $(CXXFLAGS) -I$(IDIR) -I$(EIGDIR) source/edge.cpp -o $@

build/obj/element.o : source/element.cpp
	mkdir -p build/obj
	g++ -c $(CXXFLAGS) -I$(IDIR) -I$(EIGDIR) source/element.cpp -o $@

build/obj/functions1D.o : source/functions1D.cpp
	mkdir -p build/obj
	g++ -c $(CXXFLAGS) -I$(IDIR) -I$(EIGDIR) source/functions1D.cpp -o $@

clean : 
	rm build/obj/*.o build/bin/*
