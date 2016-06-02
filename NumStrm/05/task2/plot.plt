a = 2000
#stats "res/out.dat";
set xrange[1500:1999]
set yrange[0:6000]
STATS_records = a*120
print STATS_records
do for[ii=0:STATS_records:a]{
  plot "res/out.dat" every ::ii::ii+a-1 u 1:2 w lines
  pause 0.0001
}
