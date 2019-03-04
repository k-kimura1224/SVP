#/bin/sh
#chmod u+x gen.sh

run()
{
   file=$1

   buf=${1##*/}
   name=${buf%.*}

   testval=`/usr/bin/grep "best value: " $file | sed -e 's/best value: //g'`
   bestval=`cat ./objfiles/${name}.best`

   if [ $testval = $bestval ]; then
      printf "${name}\t--->  ok!\n"
   else
      printf "${name}\t--->  !!!\n"
   fi
}

tail result.txt
echo ""

for file in `ls ./*.log`
do
   run $file
done

