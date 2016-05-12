set terminal png
set output "res/out_lf3.png"
plot "res/out_lf2.dat" u 1:3

set terminal png
set output "res/out_rk.png"
plot "res/out_RK.dat" u 1:3

set terminal png
set output "res/out_LF_old.png"
plot "res/out_LF_old.dat" u 1:3
