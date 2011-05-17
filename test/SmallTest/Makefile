all: test
PINDEL=../../src/pindel
DST=ActualOutput

test:
	-/bin/rm -r ActualOutput
	mkdir ActualOutput
	$(PINDEL) -i sim1chrVs2.conf -f sim1chrVs2.fa -c ALL -o $(DST)/sim1chrVs20305
	$(PINDEL) -p COLO-829_20-p_ok.txt -f hs_ref_chr20.fa -c 20:0-1,500,000 -o $(DST)/colowobd15
	$(PINDEL) -p COLO-829_20-p_ok.txt -f hs_ref_chr20.fa -b COLO-829.20.BreakDancer.sv -c 20:0-1,500,000  -o $(DST)/colousingbd15
	./compare.sh

clean:
	-/bin/rm -r ActualOutput
