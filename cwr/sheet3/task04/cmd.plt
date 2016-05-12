set terminal png
set output "res/out.png"
set yrange[-5:2]
plot "res/out.dat" u 1:2
