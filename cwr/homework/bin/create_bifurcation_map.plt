set ylabel "{/Symbol Q} [rad]"
set xlabel "Antriebskraft f"

set yrange[0:2.7]

set key top left

set term pngcairo enhanced
set output "../results/bifurcation/bifurcation.png"

plot "../results/bifurcation/data.dat" u 1:2  pt "." t "Ermittelte periodische Punkte in Theta"
