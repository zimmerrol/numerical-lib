set terminal pngcairo enhanced
set output "res/phase.png"
plot "res/euler.dat" u 1:2 t "Euler", "res/ab2.dat" u 1:2 t "AB2", "res/ab3.dat" u 1:2 t "AB3"
