./pendulum ../results/phase_diagrams/data_f_0.0.dat 0.3333 0.75 0.0 500 5000 1 1
./pendulum ../results/phase_diagrams/data_f_0.9.dat 0.3333 0.75 0.9 500 5000 1 1
./pendulum ../results/phase_diagrams/data_f_1.35.dat 0.3333 0.75 1.35 500 5000 1 1
./pendulum ../results/phase_diagrams/data_f_1.45.dat 0.3333 0.75 1.45 500 5000 1 1
gnuplot "create_phase_diagrams.plt"
