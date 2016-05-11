set terminal png
set output "out.png"
set yrange[-5:2]
plot "out.dat" u 1:2
