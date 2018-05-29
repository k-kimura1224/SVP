#set terminal aqua
set terminal png
set xrange[1:45]
#set yrange[0:1000]
#set key left top
set xlabel "Depth"
set ylabel "Total time"

set output "40.png"
plot "m0/BKZ20-SVP_n40_seed0.time" using 1:2 title "m0" with lines,\
"m1/BKZ20-SVP_n40_seed0.time" using 1:2 title "m1" with lines,\
"m2/BKZ20-SVP_n40_seed0.time" using 1:2 title "m2" with lines

set output "41.png"
plot "m0/BKZ20-SVP_n41_seed31.time" using 1:2 title "m0" with lines,\
"m1/BKZ20-SVP_n41_seed31.time" using 1:2 title "m1" with lines,\
"m2/BKZ20-SVP_n41_seed31.time" using 1:2 title "m2" with lines

set output "42.png"
plot "m0/BKZ20-SVP_n42_seed47.time" using 1:2 title "m0" with lines,\
"m1/BKZ20-SVP_n42_seed47.time" using 1:2 title "m1" with lines,\
"m2/BKZ20-SVP_n42_seed47.time" using 1:2 title "m2" with lines

set output "43.png"
plot "m0/BKZ20-SVP_n43_seed2.time" using 1:2 title "m0" with lines,\
"m1/BKZ20-SVP_n43_seed2.time" using 1:2 title "m1" with lines,\
"m2/BKZ20-SVP_n43_seed2.time" using 1:2 title "m2" with lines

set output "44.png"
plot "m0/BKZ20-SVP_n44_seed8.time" using 1:2 title "m0" with lines,\
"m1/BKZ20-SVP_n44_seed8.time" using 1:2 title "m1" with lines,\
"m2/BKZ20-SVP_n44_seed8.time" using 1:2 title "m2" with lines

set output "45.png"
plot "m0/BKZ20-SVP_n45_seed79.time" using 1:2 title "m0" with lines,\
"m1/BKZ20-SVP_n45_seed79.time" using 1:2 title "m1" with lines,\
"m2/BKZ20-SVP_n45_seed79.time" using 1:2 title "m2" with lines

set output "46.png"
plot "m0/BKZ20-SVP_n46_seed16.time" using 1:2 title "m0" with lines,\
"m1/BKZ20-SVP_n46_seed16.time" using 1:2 title "m1" with lines,\
"m2/BKZ20-SVP_n46_seed16.time" using 1:2 title "m2" with lines

set output "47.png"
plot "m0/BKZ20-SVP_n47_seed95.time" using 1:2 title "m0" with lines,\
"m1/BKZ20-SVP_n47_seed95.time" using 1:2 title "m1" with lines,\
"m2/BKZ20-SVP_n47_seed95.time" using 1:2 title "m2" with lines

set output "48.png"
plot "m0/BKZ20-SVP_n48_seed7.time" using 1:2 title "m0" with lines,\
"m1/BKZ20-SVP_n48_seed7.time" using 1:2 title "m1" with lines,\
"m2/BKZ20-SVP_n48_seed7.time" using 1:2 title "m2" with lines

set output "49.png"
plot "m0/BKZ20-SVP_n49_seed7.time" using 1:2 title "m0" with lines,\
"m1/BKZ20-SVP_n49_seed7.time" using 1:2 title "m1" with lines,\
"m2/BKZ20-SVP_n49_seed7.time" using 1:2 title "m2" with lines

