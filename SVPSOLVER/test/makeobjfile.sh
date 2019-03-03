#/bin/sh
#chmod u+x gen.sh

run()
{
   file=$1

   buf=${1##*/}
   name=${buf%.*}

   /usr/bin/grep "best value: " $file > ./objfiles/${name}.best
   gsed -i 's/best value: //g' ./objfiles/${name}.best
}

for file in `ls ./*.log`
do
   run $file
done

