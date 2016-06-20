#set term gif animate
#set output "animate.gif"
a = 500
stats "res/out.dat";
set xrange[0:1000000]
set yrange[-0.1:2]
print STATS_records

do for[ii=0:STATS_records:a]{
  plot "res/out.dat" every ::ii::ii+a-1 u 1:2 w lines
  #pause 0.0001
}
