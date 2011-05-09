all:
	make -C src SAMTOOLS=$(realpath $(SAMTOOLS)) all

clean:
	make -C src SAMTOOLS=$(realpath $(SAMTOOLS)) clean

test:
	make -C src SAMTOOLS=$(realpath $(SAMTOOLS)) test
