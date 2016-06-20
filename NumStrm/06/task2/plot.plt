set term png

set output "res/a1.png"
plot "res/out.dat" u 1:2 t "FCTS 1", "res/out.dat" u 1:3 t "Analytic Solution"
