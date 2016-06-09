./task01 res/out01.dat 20 20 1e-3 1 1
./task01 res/out02.dat 20 20 2e-3 1 1

./task01 res/out03.dat 20 20 1e-2 1 0.1
./task01 res/out04.dat 20 20 2e-2 1 0.1

./task01 res/out05.dat 20 20 1e-2 1 0.01

gnuplot "plot.plt"

echo "Successfully calculated and plotted all scenarios."
