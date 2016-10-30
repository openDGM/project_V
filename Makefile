IDIR = include
EIGDIR = ../eigen
CXXFLAGS = -std=c++14

all : build/bin/openDGM

build/bin/openDGM : build/obj/functions2D.o build/obj/functions1D.o build/obj/edge.o build/obj/element.o build/obj/element2D.o build/obj/edge2D.o
	mkdir -p build/bin
	g++ $(CXXFLAGS) -I$(IDIR) -I$(EIGDIR) build/obj/functions1D.o build/obj/edge.o build/obj/element.o main.cpp -o $@

build/obj/edge.o : source/edge.cpp
	mkdir -p build/obj
	g++ -c $(CXXFLAGS) -I$(IDIR) -I$(EIGDIR) source/edge.cpp -o $@

build/obj/edge2D.o : source/edge2D.cpp
	mkdir -p build/obj
	g++ -c $(CXXFLAGS) -I$(IDIR) -I$(EIGDIR) source/edge2D.cpp -o $@

build/obj/element.o : source/element.cpp
	mkdir -p build/obj
	g++ -c $(CXXFLAGS) -I$(IDIR) -I$(EIGDIR) source/element.cpp -o $@

build/obj/element2D.o : source/element2D.cpp
	mkdir -p build/obj
	g++ -c $(CXXFLAGS) -I$(IDIR) -I$(EIGDIR) source/element2D.cpp -o $@

build/obj/functions1D.o : source/functions1D.cpp
	mkdir -p build/obj
	g++ -c $(CXXFLAGS) -I$(IDIR) -I$(EIGDIR) source/functions1D.cpp -o $@

build/obj/functions2D.o : source/functions2D.cpp
	mkdir -p build/obj
	g++ -c $(CXXFLAGS) -I$(IDIR) -I$(EIGDIR) source/functions2D.cpp -o $@

clean : 
	rm build/obj/*.o build/bin/*
