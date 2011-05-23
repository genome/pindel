# Makefile for the test suites
# Include the local configuration
include ../Makefile.local

all: functional-tests coverage-tests acceptance-tests
test: functional-tests
	
clean:
	make -C SmallTest clean

acceptance-tests:
	make -C SmallTest acceptance-tests

coverage-tests:
	make -C SmallTest coverage-tests

functional-tests:
	make -C SmallTest functional-tests

regression-tests:
	echo 'No tests to run.'

# Pseudo targets for configuration
.PHONY: clean acceptance-tests coverage-tests functional-tests regression-tests
.IGNORE: clean