# Make file for pindel
#
# Enable implicit rules for .cpp
.SUFFIXES:.cpp
	
# Define the location of the samtools, you can override on the commandline
# e.g. make SAMTOOLS=../samtools-0.1.9
SAMTOOLS=samtools-0.1.9

# Define the flags for the g++ compile
# -O3 use all optimizations (including inlining of functions and unrolling loops)
# -fopenmp use openmp
# -I $(SAMTOOLS) add the samtools to the include path
# -L $(SAMTOOLS) add the samtools to the linker path 
# -lm -lz -lbam link statically with libm.a libz.a and libbam.a
CXXFLAGS=-O3 -fopenmp -I $(SAMTOOLS) -L $(SAMTOOLS) -lm -lz -lbam

all:	pindel
	
pindel:	pindel.cpp

clean:
	rm -f pindel

test:	pindel
	echo "Not implemented yet"