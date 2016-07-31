reset
set term png

set palette maxcolors 128
set xlabel "x"
set ylabel "y"
set key top right box opaque
set xrange[0:100]
set yrange[0:100]

set output "100x100real.png"
plot "out2.dat" matrix w image

set output "100x100.png"
plot "out.dat" matrix w image