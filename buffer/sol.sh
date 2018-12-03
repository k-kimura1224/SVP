#!/bin/sh

f(){
	file=$1
	buf=${1##*/}
	name=${buf%.*}

	rm -f buf.txt

	sed 's/\]//g' $file | sed 's/\[//g'  > buf.txt

	line=`cat buf.txt`

	sum=0

	for var in $line
	do
		sum=`echo "${sum} + (${var} * ${var})" | bc`
	done

	sum=`echo "scale=3; sqrt(${sum})" | bc`

	echo "${file}: $sum" >> objval.txt

	rm -f buf.txt
}

rm -f objval.txt

for file in `ls *.sol`
do
	f $file
done

cat objval.txt
