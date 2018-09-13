#/bin/sh
#chmod u+x gen.sh

run()
{
   file=$1

   buf=${1##*/}
   name=${buf%.*}

   echo ${name}
   grep "Total" $file
   grep "gap =" $file
   grep "Objective =" $file
   grep "Nodes =" $file
   echo ""
}

rm -f result.txt

for file in `ls BKZ*.log`
do
   run $file >> result.txt
done

cat result.txt

