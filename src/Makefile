# Make file for pindel
#
# Enable implicit rules for .cpp
.SUFFIXES:.cpp .xml
	
# Define the location of the samtools, you can override on the commandline
# e.g. make SAMTOOLS=../samtools-0.1.9
SAMTOOLS=/home/leon/Downloads/samtools-0.1.9

# Define the flags for the g++ compile
# -O3 use all optimizations (including inlining of functions and unrolling loops)
# -fopenmp use openmp
# -I $(SAMTOOLS) add the samtools to the include path
# -L $(SAMTOOLS) add the samtools to the linker path 
# -lm -lz -lbam link statically with libm.a libz.a and libbam.a
CXXFLAGS=-O3 -fopenmp -I$(realpath $(SAMTOOLS)) -lm -lz -L$(realpath $(SAMTOOLS)) -lbam

all:	pindel
	
pindel:	pindel.cpp
	$(CXX) $(CXXFLAGS) -O3 -fopenmp -I$(SAMTOOLS) $< -o $@ -lm -lz -L$(SAMTOOLS) -lbam

clean:
	rm -f pindel
	rm -f cppcheck-result-pindel.xml

test:	pindel cppcheck

cppcheck: cppcheck-result-pindel.xml
	
cppcheck-result-pindel.xml: pindel.cpp
	cppcheck --quiet --enable=all --xml $< 2> $@
