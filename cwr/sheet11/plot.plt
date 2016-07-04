set term gif animate
set output "animate.gif"

a = 1025
stats "out/psi_alpha_10.dat";
set xrange[-256:0]
set yrange[-0.001:.1]
print STATS_records

do for[ii=0:STATS_records-1:a]{
  plot "out/psi_alpha_10.dat" every ::ii::ii+a-1 u 1:2 w lines t sprintf("Probability alpha = 10 for t = %3d", ii/a), \
"out/psi_alpha_20.dat" every ::ii::ii+a-1 u 1:2 w lines t sprintf("Probability alpha = 20 for t = %3d", ii/a), \
"out/psi_alpha_30.dat" every ::ii::ii+a-1 u 1:2 w lines t sprintf("Probability alpha = 30 for t = %3d", ii/a), \
"out/psi_alpha_40.dat" every ::ii::ii+a-1 u 1:2 w lines t sprintf("Probability alpha = 40 for t = %3d", ii/a), \
"out/psi_alpha_50.dat" every ::ii::ii+a-1 u 1:2 w lines t sprintf("Probability alpha = 50 for t = %3d", ii/a), \
"out/psi_alpha_60.dat" every ::ii::ii+a-1 u 1:2 w lines t sprintf("Probability alpha = 60 for t = %3d", ii/a), \
"out/psi_alpha_70.dat" every ::ii::ii+a-1 u 1:2 w lines t sprintf("Probability alpha = 70 for t = %3d", ii/a), \
"out/psi_alpha_80.dat" every ::ii::ii+a-1 u 1:2 w lines t sprintf("Probability alpha = 80 for t = %3d", ii/a), \
"out/psi_alpha_90.dat" every ::ii::ii+a-1 u 1:2 w lines t sprintf("Probability alpha = 90 for t = %3d", ii/a), \
"out/psi_alpha_100.dat" every ::ii::ii+a-1 u 1:2 w lines t sprintf("Probability alpha = 100 for t = %3d", ii/a)
  #pause 0.0001
}
