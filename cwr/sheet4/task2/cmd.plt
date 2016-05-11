set terminal png
set output "out_lf3.png"
plot "out_lf2.dat" u 1:3

set terminal png
set output "out_rk.png"
plot "out_RK.dat" u 1:3

set terminal png
set output "out_LF_old.png"
plot "out_LF_old.dat" u 1:3
