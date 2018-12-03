#!/bin/sh

f(){
	file=$1
	buf=${1##*/}
	name=${buf%.*}

	echo "fplll -a svp $file"
	echo "$file" >> time.txt

	(time fplll -a svp $file > ${name}.sol) >> time.txt 2>&1

	echo  "" >> time.txt

	mv ${name}.sol solfiles
}

rm -f time.txt
rm -rf solfiles

mkdir solfiles

for file in `ls *.txt`
do
	f $file
done

cat time.txt
mv time.txt solfiles
