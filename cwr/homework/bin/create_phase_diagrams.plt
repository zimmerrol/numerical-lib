set ylabel "{/Symbol Q} [rad]"
set xlabel "{/Symbol w} [rad/s]"

set yrange[-pi:pi]
set xrange[-pi:pi]

set key top right

set term pngcairo enhanced

set output "../results/phase_diagrams/data_f_0.0.png"
plot "../results/phase_diagrams/data_f_0.0.dat" u 2:3 w lines t "Ermittelte Werte für f=0.0"

set output "../results/phase_diagrams/data_f_0.9.png"
plot "../results/phase_diagrams/data_f_0.9.dat" u 2:3 w lines t "Ermittelte Werte für f=0.9"

set output "../results/phase_diagrams/data_f_1.35.png"
plot "../results/phase_diagrams/data_f_1.35.dat" u 2:3 w lines t "Ermittelte Werte für f=1.35"

set key bottom right
set output "../results/phase_diagrams/data_f_1.45.png"
plot "../results/phase_diagrams/data_f_1.45.dat" u 2:3 w lines t "Ermittelte Werte für f=1.45"

set key bottom right
set output "../results/phase_diagrams/data_f_0.5.png"
plot "../results/phase_diagrams/data_f_0.5.dat" u 2:3 w lines t "Ermittelte Werte für f=1.5"

set key bottom right
set output "../results/phase_diagrams/data_f_0.7.png"
plot "../results/phase_diagrams/data_f_0.7.dat" u 2:3 w lines t "Ermittelte Werte für f=0.7"
