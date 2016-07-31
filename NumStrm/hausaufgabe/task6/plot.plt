set term pngcairo enhanced
set output "sor_comp.png"

set xlabel "Relaxationsparameter {/Symbol a}"
set ylabel "Iterationen"

plot "comp.dat" w lines lt 2 notitle, "comp.dat" lt 2 t "Ermittelte Iterationen"
