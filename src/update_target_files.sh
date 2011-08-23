#!/bin/bash

# place this script in the directory that contains the pindel directory

PINDEL_DIR=".."
TEST_DIR="$PINDEL_DIR/test/SmallTest"
RESULT_DIR="$TEST_DIR/TargetOutput"

$PINDEL_DIR/pindel -p $TEST_DIR/COLO-829_20-p_ok.txt \
                   -f $TEST_DIR/hs_ref_chr20.fa \
                   -b $TEST_DIR/COLO-829.20.BreakDancer.sv \
		   -c 20:0-1,500,000 -T 1\
                   -o $RESULT_DIR/colousingbd15 > coloplus.log 

$PINDEL_DIR/pindel -p $TEST_DIR/COLO-829_20-p_ok.txt \
                   -f $TEST_DIR/hs_ref_chr20.fa \
		   -c 20:0-1,500,000 -T 1\
                   -o $RESULT_DIR/colowobd15 > colominus.log

cd $TEST_DIR
../pindel -i sim1chrVs2.conf \
		   -f sim1chrVs2.fa \
                   -c ALL -T 1\
                   -o TargetOutput/sim1chrVs20305 >  ../../../../../small.log

cd -
cd pindel
svn commit --message="Updating results after pindel improvement"



