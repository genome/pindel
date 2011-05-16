# Make file for pindel
#
# Enable implicit rules for .cpp
.SUFFIXES:.cpp .o .xml
	
# Define the location of the samtools, you can override on the commandline
# e.g. make SAMTOOLS=../samtools-0.1.16
SAMTOOLS=../samtools-0.1.16
CPPCHECKS=cppcheck-result-pindel.xml cppcheck-result-reader.xml \
	  cppcheck-result-searcher.xml cppcheck-result-reporter.xml \
	  cppcheck-result-parameter.xml
OBJECTS=pindel.o reader.o reporter.o searcher.o parameter.o

all:	pindel
	
.cpp.o:
	$(CXX) $(CXXFLAGS) -Wall -c -O3 -fopenmp -I$(SAMTOOLS) $<
	
cppcheck-result-%.xml: %.cpp
	cppcheck --quiet --enable=all --xml -I$(SAMTOOLS) $< 2> $@
	
reader.o: reader.cpp pindel.h
reporter.o: reporter.cpp pindel.h
searcher.o: searcher.cpp pindel.h
pindel.o: pindel.cpp pindel.h
parameter.o: parameter.cpp parameter.h

pindel:	$(OBJECTS)
	$(CXX) $^ -O3 -fopenmp -lm -lz -L$(SAMTOOLS) -lbam -o $@

clean:
	rm -f pindel
	rm -f $(OBJECTS)
	rm -f $(CPPCHECKS)

test:	pindel cppcheck

cppcheck: $(CPPCHECKS)
