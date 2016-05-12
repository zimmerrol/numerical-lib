set terminal pngcairo enhanced
set output "res/plot.png"

plot "res/dat_out.dat" u 1:2 w lines t "Exact diff", "res/dat_out.dat" u 1:3 w lines  t "ST", "res/dat_out.dat" u 1:4 w lines  t "FW"
