reset
set term png
set palette maxcolors 128
set cbrange[0:1]
set xlabel "x"
set ylabel "y"
set key top right box opaque
set xrange[-0.5:31.5]
set yrange[-0.5:31.5]

set output "grid05.png"
plot "out05.dat" matrix w image

set output "grid005.png"
plot "out005.dat" matrix w image

set output "grid0005.png"
plot "out0005.dat" matrix w image
