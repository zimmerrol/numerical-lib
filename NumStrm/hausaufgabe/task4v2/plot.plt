reset
set term pngcairo enhanced

set palette maxcolors 128
set xlabel "x"
set ylabel "y"
set cbrange[0:1]
set key outside;
set key top center box opaque
set xrange[0:30]
set yrange[0:30]


set output "test.png"
plot "out3.dat" matrix w image
