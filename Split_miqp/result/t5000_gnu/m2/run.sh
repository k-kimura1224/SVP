#/bin/sh
#chmod u+x gen.sh

run()
{
   file=$1
   n=$(echo $(cat $file| wc -l))
   for i in `seq $n`
   do
      sed "${i}s/^/${i} /g" $file > buf.time
      cp -f buf.time $file
      rm -f buf.time
   done
}


for file in `ls ./*.time`
do
   run $file
done

