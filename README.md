# Pindel

## Compiling 

To compile Pindel you need GNU Make, GCC and cppcheck and python. Usually they 
are already installed on Linux. On the Mac you will need to install XCode (3 or
higher), the latest version can be downloaded from the Mac App Store. cppcheck 
can be installed with Fink, MacPorts or Homebrew, where homebrew is more
convenient if you don't need to install other software from source.

Pindel requires samtools; you can download the latest version of samtools from
[the Samtools Sourceforge site](http://sourceforge.net/projects/samtools/files/). 
You need to extract and build samtools before you compile Pindel. Refer to the 
documentation of samtools for the requirements to build samtools.

If you make the Pindel project the first time, it will create a `Makefile.local`.
It will try to be smart and generate this file with the correct location of
samtools if you specify it on the make commandline.

    make SAMTOOLS=~/samtools-0.1.18

The `Makefile.local` also includes a couple of options for the different tests.
The `*_TIME` options define the max execution time for the jobs in the acceptance
tests, if they are exceeded the acceptance tests will fail.

The `THREADS` option defines the number of threads used in the functional tests,
this should be at least 2, it is recommended to match the number of cores in
your system to get the fastest execution time.

## Different tests
The following test targets are defined:

1. acceptance-tests
2. coverage-tests
3. functional-tests
4. regression-tests

