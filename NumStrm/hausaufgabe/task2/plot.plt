reset
set term png
set output "grid.png"
set palette maxcolors 128
set xlabel "x"
set ylabel "y"
set cbrange[-1:1]
set key top right box opaque
set xrange[0:10]
set yrange[0:10]
plot "out.dat" matrix w image
