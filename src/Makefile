# Makefile for pindel and pindel-debug
# Include the Makefile.local from the main project for configuration
include ../Makefile.local

CPPCHECKS=cppcheck-result-pindel.xml cppcheck-result-reader.xml \
	  cppcheck-result-searcher.xml cppcheck-result-reporter.xml \
	  cppcheck-result-parameter.xml cppcheck-result-ControlState.xml
OBJECTS=pindel.o reader.o reporter.o searcher.o parameter.o ControlState.o
OBJECTS_DEBUG=pindel-debug.o reader-debug.o reporter-debug.o \
	      searcher-debug.o parameter-debug.o ControlState-debug.o

all:	pindel pindel-debug
test:	pindel $(CPPCHECKS)

# Dependencies for the normal .o files	
reader.o: reader.cpp reader.h pindel.h
reporter.o: reporter.cpp reporter.h pindel.h
searcher.o: searcher.cpp searcher.h pindel.h
pindel.o: pindel.cpp pindel.h \
	reader.h reporter.h searcher.h parameter.h ControlState.h
parameter.o: parameter.cpp parameter.h
ControlState.o: ControlState.cpp ControlState.h

# Dependencies for the debug version of the .o files
reader-debug.o: reader.cpp reader.h pindel.h
reporter-debug.o: reporter.cpp reporter.h pindel.h
searcher-debug.o: searcher.cpp searcher.h pindel.h
pindel-debug.o: pindel.cpp pindel.h \
	reader.h reporter.h searcher.h parameter.h ControlState.h
parameter-debug.o: parameter.cpp parameter.h
ControlState-debug.o: ControlState.cpp ControlState.h

# Debug version of pindel
pindel-debug: $(OBJECTS_DEBUG)
	$(CXX) $^ -fopenmp -g -pg -fprofile-arcs -ftest-coverage \
	-lm -lz -L$(SAMTOOLS) -lbam -o $@

# Normal version of pindel	
pindel:	$(OBJECTS)
	$(CXX) $^ -O3 -fopenmp -lm -lz -L$(SAMTOOLS) -lbam -o $@

clean:
	-rm -f pindel pindel-debug
	-rm -f $(OBJECTS) $(OBJECTS_DEBUG)
	-rm -f $(CPPCHECKS)
	-rm -f *.gcov *.gcda *.gcno gmon.out

# Pattern rule to create normal .o files	
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -Wall -g -c -O3 -fopenmp -I$(SAMTOOLS) $< -o $@

# Pattern rule to create debug version of .o files
%-debug.o: %.cpp
	$(CXX) $(CXXFLAGS) -Wall -g -c -fopenmp -g -pg \
	-fprofile-arcs -ftest-coverage -I$(SAMTOOLS) $< -o $@
	
# Pattern rule to create cppcheck files	
cppcheck-result-%.xml: %.cpp
	cppcheck --quiet --enable=all --xml -I$(SAMTOOLS) $(realpath $<) 2> $@

# Pseudo targets for configuration
.PHONY: all clean test
.IGNORE: clean