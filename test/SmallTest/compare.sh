#!/bin/sh
for f in sim1chrVs20305_{BP,D,INV,LI,SI,TD}
do
	diff -c TargetOutput/$f ActualOutput/$f || (echo Not identical ; exit 1)
done || (echo Something went wrong ; exit 1)
echo Results identical to expected output. OK.
exit 0
