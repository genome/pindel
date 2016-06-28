#!/bin/sh

cd ..
rm src/*.o
rm src/*.d
rm pindel
rm pindel2vcf
rm pindel2vcf4tcga
rm sam2pindel
rm Makefile.local

./INSTALL ../htslib 2>&1 | tee devtools/installer_output.txt

