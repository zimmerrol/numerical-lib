#set term gif animate
#set output "animate.gif"
a = 1025
stats "out/out.dat";
set xrange[-256:0]
set yrange[-0.001:.1]
print STATS_records

do for[ii=0:STATS_records-1:a]{
  plot "out/out.dat" every ::ii::ii+a-1 u 1:2 w lines
  #pause 0.0001
}
