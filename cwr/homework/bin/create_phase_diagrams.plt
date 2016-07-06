set ylabel "{/Symbol Q} [rad]"
set xlabel "{/Symbol w} [rad/s]"

set yrange[-pi:pi]
set xrange[-pi:pi]

set key bottom right

set term pngcairo enhanced

set output "../results/phase_diagrams/data_f_1.35.png"
plot "../results/phase_diagrams/data_f_1.35.dat" u 2:3 w lines t "Ermittelte Werte für f=1.35"

set key bottom right
set output "../results/phase_diagrams/data_f_1.40.png"
plot "../results/phase_diagrams/data_f_1.40.dat" u 2:3 w lines t "Ermittelte Werte für f=1.40"
