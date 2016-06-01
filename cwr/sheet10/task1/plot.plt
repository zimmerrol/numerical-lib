a = 100
stats "res/out.dat";
set xrange[0:100]
set yrange[-3:3]
print STATS_records
do for[ii=0:STATS_records:a]{
  plot "res/out.dat" every ::ii::ii+a-1 u 1:2 w lines
  pause 0.0001
}
