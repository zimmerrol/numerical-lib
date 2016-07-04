set ylabel "{/Symbol Q} [rad]"
set y2label "{/Symbol w} [rad/s]"
set xlabel "Zeit t"

set key top right


set y2tics
set autoscale  y
set autoscale  y2

set term pngcairo enhanced

set output "../results/theta_omega_tn/data_f_0.0.png"
plot "../results/theta_omega_tn/data_f_0.0.dat" u 1:2 w lines t "Omega" axes x1y1, "../results/theta_omega_tn/data_f_0.0.dat" u 1:3 w lines t "Theta" axes x1y2

set output "../results/theta_omega_tn/data_f_0.9.png"
plot "../results/theta_omega_tn/data_f_0.9.dat" u 1:2 w lines t "Omega" axes x1y1, "../results/theta_omega_tn/data_f_0.9.dat" u 1:3 w lines t "Theta" axes x1y2

set output "../results/theta_omega_tn/data_f_1.35.png"
plot "../results/theta_omega_tn/data_f_1.35.dat" u 1:2 w lines t "Omega" axes x1y1, "../results/theta_omega_tn/data_f_1.35.dat" u 1:3 w lines t "Theta" axes x1y2

set key bottom right
set output "../results/theta_omega_tn/data_f_1.45.png"
plot "../results/theta_omega_tn/data_f_1.45.dat" u 1:2 w lines t "Omega" axes x1y1, "../results/theta_omega_tn/data_f_1.45.dat" u 1:3 w lines t "Theta" axes x1y2
