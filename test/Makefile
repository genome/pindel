all: functional-tests coverage-tests acceptance-tests
test: functional-tests
	
.PHONY : test clean

clean:
	make -C SmallTest clean

acceptance-tests:
	make -C SmallTest acceptance-tests

coverage-tests:
	make -C SmallTest coverage-tests

functional-tests:
	make -C SmallTest functional-tests
	