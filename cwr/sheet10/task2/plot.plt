a = 100
#stats "res/out.dat";
STATS_records = 1000*100
set xrange[0:500]
set yrange[-5:10]
print STATS_records
do for[ii=0:STATS_records:a]{
  print "aa"
  plot "res/out.dat" every ::ii::ii+a-1 u 1:2 w lines
  print "bb"
  #pause 0.0001
}
