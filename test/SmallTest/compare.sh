#!/bin/bash

function exitWithMessage() {
  echo "$1"
  exit 1
}

for t in colowobd15 sim1chrVs20305 colousingbd15
do
  for f in ${t}_{BP,D,INV,LI,SI,TD}
  do
	  diff -c TargetOutput/$f ActualOutput/$f || exitWithMessage 'Not identical'
  done || exitWithMessage 'Inner for loop failed'
done || exitWithMessage 'Outer for loop failed'
echo Results identical to expected output. OK.
exit 0
