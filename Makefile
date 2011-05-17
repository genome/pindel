# Define the location of the samtools, you can override on the commandline
# e.g. make SAMTOOLS=~/samtools-0.1.16
SAMTOOLS=samtools-0.1.16
.PHONY : test all test coverage

all:  test coverage

pindel:
	make -C src SAMTOOLS=$(realpath $(SAMTOOLS)) pindel
clean:
	make -C src SAMTOOLS=$(realpath $(SAMTOOLS)) clean
	make -C test clean

test:
	make -C src SAMTOOLS=$(realpath $(SAMTOOLS)) test
	make -C test test

coverage:
	make -C src SAMTOOLS=$(realpath $(SAMTOOLS)) pindel-debug
	make -C test coverage