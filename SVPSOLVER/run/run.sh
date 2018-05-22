#/bin/sh
#chmod u+x gen.sh

run()
{
   file=$1

   buf=${1##*/}
   name=${buf%.*}

   thread=$2
   time=$3

   echo "../bin/SVPSOLVER $file $thread $time > ${name}.log"
   ../bin/SVPSOLVER $file $thread $time > ${name}.log
}

dirname=$1
thread=$2
time=$3

for file in `ls ${dirname}/*.dat`
do
   run $file $thread $3
done

