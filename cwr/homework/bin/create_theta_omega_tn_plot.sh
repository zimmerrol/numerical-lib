./pendulum ../results/theta_omega_tn/data_f_0.0.dat 0.3333 0.75 0.0 500 5000
./pendulum ../results/theta_omega_tn/data_f_0.9.dat 0.3333 0.75 0.9 500 5000
./pendulum ../results/theta_omega_tn/data_f_1.35.dat 0.3333 0.75 1.35 500 5000
./pendulum ../results/theta_omega_tn/data_f_1.45.dat 0.3333 0.75 1.45 500 5000
gnuplot "create_theta_omega_tn_plot.plt"
