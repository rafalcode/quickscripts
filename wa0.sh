#!/bin/bash
EXPECTED_ARGS=2 # change value to suit!
if [ $# -ne $EXPECTED_ARGS ]; then
	echo "Correct usage: $0 <args ..>"
	exit
fi
T=$1
Q=$2
OF=al_${T%.*}_${Q%.*}.fa
GO=10.0
GE=2.5
echo "water -asequence $T -bsequence $Q -gapopen $GO -gapextend $GE -outfile $OF -aformat3 fasta"
water -asequence $T -bsequence $Q -gapopen $GO -gapextend $GE -outfile $OF -aformat3 fasta
