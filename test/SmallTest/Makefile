all: test

test:
	-/bin/rm -r ActualOutput
	mkdir ActualOutput
	../../src/pindel -i sim1chrVs2.conf -f sim1chrVs2.fa -c ALL -o ActualOutput/sim1chrVs20305
	diff -cr TargetOutput ActualOutput || (echo Not identical ; exit 1)

clean:
	-/bin/rm -r ActualOutput
