# Make file for pindel
#
# Enable implicit rules for .cpp
.SUFFIXES:.cpp .o .xml
	
# Define the location of the samtools, you can override on the commandline
# e.g. make SAMTOOLS=../samtools-0.1.16
SAMTOOLS=../samtools-0.1.16

all:	pindel
	
.cpp.o:
	$(CXX) $(CXXFLAGS) -c -O3 -fopenmp -I$(SAMTOOLS) $<

reader.o: reader.cpp pindel.h
reporter.o: reporter.cpp pindel.h
searcher.o: searcher.cpp pindel.h
pindel.o: pindel.cpp pindel.h

pindel:	pindel.o reader.o reporter.o searcher.o
	$(CXX) $^ -O3 -fopenmp -lm -lz -L$(SAMTOOLS) -lbam -o $@

clean:
	rm -f pindel reader.o reporter.o searcher.o pindel.o
	rm -f cppcheck-result-pindel.xml

test:	pindel cppcheck

cppcheck: cppcheck-result-pindel.xml
	
cppcheck-result-pindel.xml: pindel.cpp
	cppcheck --quiet --enable=all --xml $< 2> $@
