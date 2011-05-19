# Make file for pindel
#
.PHONY: all clean test
	
# Define the location of the samtools, you can override on the commandline
# e.g. make SAMTOOLS=../samtools-0.1.16
SAMTOOLS=../samtools-0.1.16
CPPCHECKS=cppcheck-result-pindel.xml cppcheck-result-reader.xml \
	  cppcheck-result-searcher.xml cppcheck-result-reporter.xml \
	  cppcheck-result-parameter.xml cppcheck-result-ControlState.xml
OBJECTS=pindel.o reader.o reporter.o searcher.o parameter.o ControlState.o
OBJECTS_DEBUG=pindel-debug.o reader-debug.o reporter-debug.o \
	      searcher-debug.o parameter-debug.o ControlState-debug.o

all:	pindel pindel-debug
	
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -Wall -c -O3 -fopenmp -I$(SAMTOOLS) $< -o $@
	
cppcheck-result-%.xml: %.cpp
	cppcheck --quiet --enable=all --xml -I$(SAMTOOLS) $< 2> $@
	
reader.o: reader.cpp pindel.h
reporter.o: reporter.cpp pindel.h
searcher.o: searcher.cpp pindel.h
pindel.o: pindel.cpp pindel.h
parameter.o: parameter.cpp parameter.h
ControlState.o: ControlState.cpp ControlState.h

%-debug.o: %.cpp
	$(CXX) $(CXXFLAGS) -Wall -c -fopenmp -g -pg \
	-fprofile-arcs -ftest-coverage -I$(SAMTOOLS) $< -o $@
	
reader-debug.o: reader.cpp pindel.h
reporter-debug.o: reporter.cpp pindel.h
searcher-debug.o: searcher.cpp pindel.h
pindel-debug.o: pindel.cpp pindel.h
parameter-debug.o: parameter.cpp parameter.h
ControlState-debug.o: ControlState.cpp ControlState.h

pindel-debug: $(OBJECTS_DEBUG)
	$(CXX) $^ -fopenmp -g -pg -fprofile-arcs -ftest-coverage \
	-lm -lz -L$(SAMTOOLS) -lbam -o $@
	
pindel:	$(OBJECTS)
	$(CXX) $^ -O3 -fopenmp -lm -lz -L$(SAMTOOLS) -lbam -o $@

clean:
	-rm -f pindel pindel-debug
	-rm -f $(OBJECTS) $(OBJECTS_DEBUG)
	-rm -f $(CPPCHECKS)
	-rm -f *.gcov *.gcda *.gcno gmon.out

test:	pindel cppcheck

cppcheck: $(CPPCHECKS)
