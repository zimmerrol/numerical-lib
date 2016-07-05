set ylabel "{/Symbol Q} [rad]"
set xlabel "Antriebskraft f"

set yrange[0:2.7]

set key top left

set term pngcairo enhanced
set output "../results/bifurcation/bifurcation.png"

plot "../results/bifurcation/data.dat" u 1:($3 < pi ? $3+0.3 : $3-2*pi+0.3)  pt "." t "Ermittelte periodische Punkte in Theta"
