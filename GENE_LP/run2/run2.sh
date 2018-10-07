#!/bin/sh

blk=20

#for i in 60 61 62 63 64 65 66 67 68 69 70
#do
#	for j in 0 1 2 3 4 5 6 7 8 9
#	do
#		$run $i $j $blk
#	done
#done

#for i in  51 52 53 54 55 56 57 58 59
#do
#	for j in 0 1 2 3 4 5 6 7 8 9
#	do
#		$run $i $j $blk
#	done
#done

for i in  40 43 46 49 
do
	for j in 5 6 7 8 9
	do
		rm -f BKZ20-SVP_n${i}_seed${j}.lp
	done
done
