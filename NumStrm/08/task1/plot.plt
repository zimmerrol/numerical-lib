set term png

set output "res/a1.png"
plot "res/out01.dat" u 1:2 t "FCTS 1", "res/out01.dat" u 1:3 t "Analytic Solution"

set output "res/a2.png"
plot "res/out02.dat" u 1:2 t "FCTS 2", "res/out02.dat" u 1:3 t "Analytic Solution"

set output "res/b1.png"
plot "res/out03.dat" u 1:2 t "FCTS 1", "res/out03.dat" u 1:3 t "Analytic Solution"

set output "res/b2.png"
plot "res/out04.dat" u 1:2 t "FCTS 2", "res/out04.dat" u 1:3 t "Analytic Solution"

set output "res/c.png"
plot "res/out05.dat" u 1:2 t "FCTS", "res/out05.dat" u 1:3 t "Analytic Solution"
