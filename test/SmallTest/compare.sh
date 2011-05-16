#!/bin/bash
for f in sim1chrVs20305_{BP,D,INV,LI,SI,TD}
do
	if ! diff -c TargetOutput/$f ActualOutput/$f
        then
            echo 'Not identical'
            exit 1
        fi
done
if test $? -ne 0
then
  echo 'For loop failed'
  exit 1
fi
echo Results identical to expected output. OK.
exit 0
