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

set output "2_05.png"
plot "2_05.dat" matrix w image t "Temperatur T zu t = 0.5"
set output "2_005.png"
plot "2_005.dat" matrix w image t "Temperatur T zu t = 0.05"
set output "2_0005.png"
plot "2_0005.dat" matrix w image t "Temperatur T zu t = 0.005"

set output "10_05.png"
plot "10_05.dat" matrix w image t "Temperatur T zu t = 0.5"
set output "10_005.png"
plot "10_005.dat" matrix w image t "Temperatur T zu t = 0.05"
set output "10_0005.png"
plot "10_0005.dat" matrix w image t "Temperatur T zu t = 0.005"
