set term png
set output "out_y.png"

set arrow from 21,0 to 21,4.66 nohead
plot "out.dat" u 2:3

reset
set term png
set output "out_z.png"

set arrow from 21,0 to 21,2.44 nohead
plot "out.dat" u 2:4

set term png
set output "out2_y.png"

set arrow from 21,0 to 21,4.66 nohead
plot "out2.dat" u 2:3

reset
set term png
set output "out2_z.png"

set arrow from 21,0 to 21,2.44 nohead
plot "out2.dat" u 2:4
