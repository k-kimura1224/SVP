#/bin/sh
#chmod u+x gen.sh

run()
{
   file=$1
   buf=${1##*/}
   name=${buf%.*}
   a='"aaa"$a'
   echo $a

   echo "set terminal png" >> t5000.plt
   echo '-n set output"' >> t5000.plt
   echo -n "${name}.png" >> t5000.plt
   echo '"' >> t5000.plt

   #rm -f t5000.plt
}

rm -f t5000.plt

for file in `ls m0/*.time`
do
   run $file
done

