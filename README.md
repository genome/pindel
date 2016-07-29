# Pindel

## Compiling 

To compile Pindel you need three things: GNU Make and GCC (which usually are 
already installed on Linux) and htslib. The last is not installed on Linux
by default, but it can be retrieved with:

git clone https://github.com/samtools/htslib

htslib needs to be built before you can start installing Pindel.
(Go to htslib's directory, and follow the directions
in its 'readme'. At the time of writing this (February 2016) it simply works
if you give the commands "make" and "sudo make install".)

To compile Pindel on OS X, you may need to do more work - 'regular' 
installation under OS X may work, but in some cases OS X gives problems with
the OpenMP library Pindel uses for speedup. In those cases, please follow the instruction on the following page to update your gcc

http://www.ficksworkshop.com/blog/14-coding/65-installing-gcc-on-mac

If htslib has been cloned and installed, go to the pindel directory 
([my-path]/pindel) and use the INSTALL script there in the following way:

./INSTALL [path-to-htslib]

for example

./INSTALL ../htslib

After this, you can run pindel by using

./pindel [options]

Plain "./pindel" without command line options will list all available command
line options, the FAQ in the Pindel root directory includes a usage example.

If there are any problems with installing or running Pindel, you may be
able to find the solution in the FAQ (the FAQ file stored in the same
directory as the INSTALL script), otherwise, feel free to open an issue 
on github (https://github.com/genome/pindel/issues) or to contact the
main author, Kai Ye, at kaiye@xjtu.edu.cn

