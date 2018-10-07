#/bin/sh
#chmod u+x gen.sh

run()
{
   file=$1

   buf=${1##*/}
   name=${buf%.*}

   time=$2

   rm -f run.txt
   echo "set timelimit ${time}" >> run.txt
   echo "set mip tolerances mipgap 1e-10" >> run.txt
   echo "read ${file}" >> run.txt
   echo "opt" >> run.txt
   echo "write ${name}.sol" >> run.txt
   echo "q" >> run.txt

   echo "Solving ${name} via CPLEX"
   ./cplex.sh > ${name}.log
}

time=$1

for file in `ls *.lp`
do
   run $file $time
done

rm -f run.txt
