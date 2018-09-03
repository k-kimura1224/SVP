#!/bin/sh
dim=$1
seed=$2
blk=$3
head=SVP_n${dim}_seed${seed}
ext=txt
./generate_random --dim ${dim} --seed ${seed} > ${head}.${ext}
echo "./generate_random --dim ${dim} --seed ${seed}"
tmpfile=tmpfile
sed -e "s/[\r\n]\+//g" ${head}.${ext} > ${tmpfile}
#tmp2=`echo $tmp | tr ' ' ', '`
#tmp2=`echo $tmp | sed -e "s/\] \[/;/g"`
#tmp3=`echo $tmp2 | sed -e "s/ /,/g"`
#tmp4=`echo $tmp3 | sed -e "s/\[\[/\[/g"`
#tmp5=`echo $tmp4 | sed -e "s/\],\]/\]/g"`
#echo $tmp
#echo $tmp2
#echo $tmp3
#echo $tmp4
dfile=BKZ${blk}-${head}.dat
if [ -e ${dfile} ]; then
	rm ${dfile}
fi
touch ${dfile}
#cat ${tmpfile}
#cat ${dfile}
fplll -a bkz -b ${blk} ${tmpfile} > ${dfile}
echo "fplll -a bkz -b ${blk}"
sed 's/\]//g' ${dfile} | sed 's/\[//g'  | sed 's/,//g' > ${tmpfile}
mv ${tmpfile} ${dfile}
echo ${dim} > ${tmpfile}
cat ${dfile} >> ${tmpfile}
mv ${tmpfile} ${dfile}
nkf -Lu ${dfile} > ${tmpfile}
mv ${tmpfile} ${dfile}
