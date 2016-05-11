set output "res/phaseRK4.png"
set terminal png
plot "out.dat" u 6:7 t "RK4"

set output "res/phaseRK2.png"
set terminal png
plot "out.dat" u 4:5 t "RK2"

set output "res/phaseEU.png"
set terminal png
plot "out.dat" u 2:3 t "EU"

set output "res/xvEU.png"
set terminal png
plot "out.dat" u 1:2 t "x", "out.dat" u 1:3 t "v"

set output "res/xvRK2.png"
set terminal png
plot "out.dat" u 1:4 t "x", "out.dat" u 1:5 t "v"

set output "res/xvRK4.png"
set terminal png
plot "out.dat" u 1:6 t "x", "out.dat" u 1:7 t "v"

set output "res/energy.png"
set terminal png
plot "out.dat" u 1:($3*$3*0.5+0.5*$2*$2) t "EU", "out.dat" u 1:($5*$5*0.5+0.5*$4*$4) t "RK2", "out.dat" u 1:($7*$7*0.5+0.5*$6*$6) t "RK4"
