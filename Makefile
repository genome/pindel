# Define the location of the samtools, you can override on the commandline
# e.g. make SAMTOOLS=~/samtools-0.1.16
SAMTOOLS=samtools-0.1.16

all:
	make -C src SAMTOOLS=$(realpath $(SAMTOOLS)) all

clean:
	make -C src SAMTOOLS=$(realpath $(SAMTOOLS)) clean

test:
	make -C src SAMTOOLS=$(realpath $(SAMTOOLS)) test
