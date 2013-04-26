# Makefile for the SmallTest testsuite
# Include the local configuration
include ../../Makefile.local

# Local configuration
PINDEL=$(realpath ../../src/pindel)
PINDEL_DBG=$(realpath ../../src/pindel-debug)
SRCDIR=$(realpath ../../src)
GCOVR=$(realpath gcovr)
TIMEX=$(realpath ./execution_time.sh)
DST=ActualOutput
COVFILES=$(SRCDIR)/pindel.cpp.gcov $(SRCDIR)/searcher.cpp.gcov \
	 $(SRCDIR)/reader.cpp.gcov $(SRCDIR)/reporter.cpp.gcov \
	 $(SRCDIR)/parameter.cpp.gcov
COVERAGE=$(abspath ./coverage.xml)
COLOUSINGBD15_ARGS=-p COLO-829_20-p_ok.txt -f hs_ref_chr20.fa -l -k -B 100 -u 0.05 -x 5\
		   -b COLO-829.20.BreakDancer.sv -Q colo_bd_out.txt -c 20:0-1,500,000 \
		   -o $(DST)/colousingbd15
COLOWOBD15_ARGS=-P pindel-config-file.txt -f hs_ref_chr20.fa -l -k -B 100 -u 0.05 -x 5 -c 20:0-1,500,000 \
	        -o $(DST)/colowobd15
SIM1CHRVS20305_ARGS=-i sim1chrVs2.conf -f sim1chrVs2.fa -c ALL \
		    -o $(DST)/sim1chrVs20305

all: test
test: functional-tests coverage-tests acceptance-tests
functional-tests: small-functional-test
acceptance-tests: execution-time
coverage-tests: $(COVERAGE)

# Functional test will verify that the results still match
small-functional-test:
	-rm -rf $(DST)
	mkdir -p $(DST)
	$(PINDEL) -T 1 $(SIM1CHRVS20305_ARGS) > sim1chrVs20305.log
	$(PINDEL) -T 1 $(COLOWOBD15_ARGS) > colowobd15.log
	$(PINDEL) -T 1 $(COLOUSINGBD15_ARGS) > colousingbd15.log
	./compare.sh

# Verify that the time to execute will not increase
execution-time:
	-rm -rf $(DST)
	mkdir -p $(DST)
	# Run in 1 thread otherwise it will skew on multicore machines
	$(TIMEX) $(SIM1CHRVS20305_TIME) $(PINDEL) -T 1 $(SIM1CHRVS20305_ARGS) > timex-sim1chrVs20305.log
	$(TIMEX) $(COLOWOBD15_TIME) $(PINDEL) -T 1 $(COLOWOBD15_ARGS) > timex-colowobd15.log
	$(TIMEX) $(COLOUSINGBD15_TIME) $(PINDEL) -T 1 $(COLOUSINGBD15_ARGS) > timex-colousingbd15.log

# Generate code coverage reports
$(COVERAGE): assemble-coverage
	cd $(SRCDIR) ; $(GCOVR) -r $(realpath $(SRCDIR)) --xml > $@

# Run the tests to get coverage	
assemble-coverage:
	-rm -rf $(DST)
	mkdir -p $(DST)
	# Run in 1 thread otherwise it will fail 
	$(PINDEL_DBG) -T 1 $(SIM1CHRVS20305_ARGS) > cover-sim1chrVs20305.log
	$(PINDEL_DBG) -T 1 $(COLOWOBD15_ARGS) > cover-colowobd15.log
	$(PINDEL_DBG) -T 1 $(COLOUSINGBD15_ARGS) > cover-colousingbd15.log
	#$(PINDEL_DBG) -h > /dev/null

clean:
	-rm -rf $(DST) $(COVDST) $(COVERAGE) gmon.out *.log *.order *.sum

gmon.sum:
	make -C ../../src CXXFLAGS=-pg clean pindel
	$(PINDEL) -T $(THREADS) $(SIM1CHRVS20305_ARGS) > cover-sim1chrVs20305.log
	gprof -b -s $(PINDEL) gmon.out
	
# Pseudo targets for configuration
.PHONY: all clean test pindel pindel-debug cppcheck acceptance-tests \
	coverage-tests functional-tests gmon.sum
.IGNORE: clean
