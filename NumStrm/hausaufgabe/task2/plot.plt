reset
set term pngcairo enhanced

set palette maxcolors 256
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

dx = 1.0/31.0
dy = 1.0/31.0
vx(x,y) = (x*dx < 0.95) ? (x*dx>0.05 ? pi*sin(2*pi*x*dx)*cos(pi*y*dy)/2 : -10) : 10
vy (x,y)= (y*dy < 0.95) ? (y*dy>0.05 ? -2*pi*cos(2*pi*x*dx)*sin(pi*y*dy)/2 : -10) : 10

set xrange[0:30]
set yrange[0:30]
set samples 15
set isosamples 15
set output "10_05_velocity.png"
set key top center horizontal nobox opaque
plot "10_05.dat" matrix w image t "Temperatur T zu t = 0.5", "++" using 1:2:(vx($1, $2)):(vy($1, $2)) w vectors filled head lt 3 lc rgb "black" t "Geschwindigkeitsvektorfeld"
#"out3.dat" matrix w image t "test",
