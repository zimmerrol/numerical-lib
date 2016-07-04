set ylabel "{/Symbol Q} [rad]"
set xlabel "Antriebskraft f"

set yrange[2.7:5.2]

set key top left

set term pngcairo enhanced
set output "../results/bifurcation/bifurcation.png"

plot "../results/bifurcation/data.dat" u 1:3 t "Ermittelte Fixpunkte in Theta"
