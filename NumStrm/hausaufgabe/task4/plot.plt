reset
set term png
set output "grid.png"
set palette maxcolors 128
set xlabel "x"
set ylabel "y"
set key top right box opaque
set xrange[-0.5:31.5]
set yrange[-0.5:31.5]
plot "out.dat" matrix w image
