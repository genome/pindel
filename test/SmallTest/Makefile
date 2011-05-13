all: test

test:
	-/bin/rm -r ActualOutput
	mkdir ActualOutput
	../../src/pindel -i sim1chrVs2.conf -f sim1chrVs2.fa -c ALL -o ActualOutput/sim1chrVs20305
	./compare.sh

clean:
	-/bin/rm -r ActualOutput
