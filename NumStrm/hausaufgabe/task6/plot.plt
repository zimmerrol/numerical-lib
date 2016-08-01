set term pngcairo enhanced


set xlabel "Relaxationsparameter {/Symbol a}"
set ylabel "Iterationen"

set output "sor_comp.png"
plot "comp.dat" w lines lt 2 notitle, "comp.dat" lt 2 t "Ermittelte Iterationen"

set output "speed_comp.png"
plot "speed_comp.dat" w lines lt 2 notitle, "speed_comp.dat" lt 2 t "Explizites Verfahren (FTCS)", "speed_comp.dat" u 1:3 w lines lt 3 notitle, "speed_comp.dat" u 1:3 lt 3 t "Implizites Verfahren (BTCS)"
