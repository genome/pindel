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
# It doesn't harm if this is larger than the actual CPUs
THREADS=4

# Targets that do not represent actual files
.PHONY: all test functional-tests acceptance-tests coverage-tests \
	execution-time assemble-coverage clean

all: test
test: functional-tests coverage-tests acceptance-tests
functional-tests: small-functional-test
acceptance-tests: execution-time
coverage-tests: $(COVERAGE)

# Functional test will verify that the results still match
small-functional-test:
	-/bin/rm -r $(DST)
	mkdir $(DST)
	$(PINDEL) -T $(THREADS) -i sim1chrVs2.conf -f sim1chrVs2.fa -c ALL -o $(DST)/sim1chrVs20305 > sim1chrVs20305.log
	$(PINDEL) -T $(THREADS) -p COLO-829_20-p_ok.txt -f hs_ref_chr20.fa -c 20:0-1,500,000 -o $(DST)/colowobd15 > colowobd15.log
	$(PINDEL) -T $(THREADS) -p COLO-829_20-p_ok.txt -f hs_ref_chr20.fa -b COLO-829.20.BreakDancer.sv -c 20:0-1,500,000  -o $(DST)/colousingbd15 > colousingbd15.log
	./compare.sh

# Verify that the time to execute will not increase
execution-time:
	-/bin/rm -r $(DST)
	mkdir $(DST)
	# Run in 1 thread otherwise it will skew on multicore machines
	$(TIMEX) 22 $(PINDEL) -T 1 -i sim1chrVs2.conf -f sim1chrVs2.fa -c ALL -o $(DST)/sim1chrVs20305 > timex-sim1chrVs20305.log
	$(TIMEX) 26 $(PINDEL) -T 1 -p COLO-829_20-p_ok.txt -f hs_ref_chr20.fa -c 20:0-1,500,000 -o $(DST)/colowobd15 > timex-colowobd15.log
	$(TIMEX) 27 $(PINDEL) -T 1 -p COLO-829_20-p_ok.txt -f hs_ref_chr20.fa -b COLO-829.20.BreakDancer.sv -c 20:0-1,500,000  -o $(DST)/colousingbd15 > timex-colousingbd15.log

# Generate code coverage reports
$(COVERAGE): assemble-coverage
	cd $(SRCDIR) ; $(GCOVR) -r . --xml > $@

# Run the tests to get coverage	
assemble-coverage:
	-/bin/rm -r $(DST)
	mkdir -p $(DST)
	# Run in 1 thread otherwise it will fail 
	$(PINDEL_DBG) -T 1 -i sim1chrVs2.conf -f sim1chrVs2.fa -c ALL -o $(DST)/sim1chrVs20305 > cover-sim1chrVs20305.log
	$(PINDEL_DBG) -T 1 -p COLO-829_20-p_ok.txt -f hs_ref_chr20.fa -c 20:0-1,500,000 -o $(DST)/colowobd15 > cover-sim1chrVs20305.log
	$(PINDEL_DBG) -T 1 -p COLO-829_20-p_ok.txt -f hs_ref_chr20.fa -b COLO-829.20.BreakDancer.sv -c 20:0-1,500,000  -o $(DST)/colousingbd15 > cover-colousingbd15.log

clean:
	-/bin/rm -r $(DST) $(COVDST) $(COVERAGE) gmon.out *.log
