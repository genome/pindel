PINDEL=../../src/pindel
PINDEL_DBG=../../src/pindel-debug
SRCDIR=$(realpath ../../src)
GCOVR=$(realpath gcovr)
DST=ActualOutput
COVFILES=$(SRCDIR)/pindel.cpp.gcov $(SRCDIR)/searcher.cpp.gcov \
	 $(SRCDIR)/reader.cpp.gcov $(SRCDIR)/reporter.cpp.gcov \
	 $(SRCDIR)/parameter.cpp.gcov
COVERAGE=$(abspath ./coverage.xml)
COVDST=tmp
.PHONY: all test coverage coverage-test assemble-coverage

all: test coverage

test:
	-/bin/rm -r $(DST)
	mkdir $(DST)
	$(PINDEL) -i sim1chrVs2.conf -f sim1chrVs2.fa -c ALL -o $(DST)/sim1chrVs20305
	$(PINDEL) -p COLO-829_20-p_ok.txt -f hs_ref_chr20.fa -c 20:0-1,500,000 -o $(DST)/colowobd15
	$(PINDEL) -p COLO-829_20-p_ok.txt -f hs_ref_chr20.fa -b COLO-829.20.BreakDancer.sv -c 20:0-1,500,000  -o $(DST)/colousingbd15
	./compare.sh

coverage: $(COVERAGE)

$(COVERAGE): coverage-test
	cd $(SRCDIR) ; $(GCOVR) -r . --xml > $@
	
coverage-test:
	-/bin/rm -r $(COVDST)
	mkdir -p $(COVDST)
	$(PINDEL_DBG) -i sim1chrVs2.conf -f sim1chrVs2.fa -c ALL -o $(COVDST)/sim1chrVs20305
	$(PINDEL_DBG) -p COLO-829_20-p_ok.txt -f hs_ref_chr20.fa -c 20:0-1,500,000 -o $(COVDST)/colowobd15
	$(PINDEL_DBG) -p COLO-829_20-p_ok.txt -f hs_ref_chr20.fa -b COLO-829.20.BreakDancer.sv -c 20:0-1,500,000  -o $(COVDST)/colousingbd15
	
clean:
	-/bin/rm -r $(DST) $(COVDST) $(COVERAGE) gmon.out
