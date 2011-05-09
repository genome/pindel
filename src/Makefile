# Make file for pindel
#
# Enable implicit rules for .cpp
.SUFFIXES:.cpp .xml
	
# Define the location of the samtools, you can override on the commandline
# e.g. make SAMTOOLS=../samtools-0.1.9
SAMTOOLS=samtools-0.1.9

# Define the flags for the g++ compile on the target itself, since it fails otherwise.
# -O3 use all optimizations (including inlining of functions and unrolling loops)
# -fopenmp use openmp
# -I $(SAMTOOLS) add the samtools to the include path
# -L $(SAMTOOLS) add the samtools to the linker path 
# -lm -lz -lbam link statically with libm.a libz.a and libbam.a
#  g++ -O3 -fopenmp -I /opt/samtools-0.1.9/ pindel.cpp -o pindel -lm -lz -L/opt/samtools-0.1.9/ -lbam
CXXFLAGS=

all:	pindel
	
pindel:	pindel.cpp
	$(CXX) -O3 -fopenmp -I $(realpath $(SAMTOOLS)) $< -o $@ -lm -lz -L$(realpath $(SAMTOOLS)) -lbam

clean:
	rm -f pindel
	rm -f cppcheck-result-pindel.xml

test:	pindel cppcheck

cppcheck: cppcheck-result-pindel.xml
	
cppcheck-result-pindel.xml: pindel.cpp
	cppcheck --quiet --enable=all --xml $< 2> $@
