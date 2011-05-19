# Define the location of the samtools, you can override on the commandline
# e.g. make SAMTOOLS=~/samtools-0.1.16
SAMTOOLS=samtools-0.1.16
.PHONY : test all test coverage

all: pindel cppcheck functional-tests coverage-tests acceptance-tests
test: pindel cppcheck functional-tests

pindel:
	make -C src SAMTOOLS=$(realpath $(SAMTOOLS)) pindel

pindel-debug:
	make -C src SAMTOOLS=$(realpath $(SAMTOOLS)) pindel-debug

cppcheck: pindel
	make -C src SAMTOOLS=$(realpath $(SAMTOOLS)) test

acceptance-tests: pindel
	make -C test acceptance-tests

coverage-tests: pindel-debug
	make -C test coverage-tests

functional-tests: pindel
	make -C test functional-tests

clean:
	make -C src SAMTOOLS=$(realpath $(SAMTOOLS)) clean
	make -C test clean
