set grid
#set size square
set xrange [-400:400]
set yrange [-400:400]


#Erwartetes Format: [Zeit]	[SatX]	[SatY]	[EarthX]	[EarthY]	[MoonX]	[MoonY]

#Hiermit wird die Zahl die Eintraege ausgelesen
stats 'res/out2.dat'

#ii laeuft ueber alle Datenpunkte in 200er-Schritten
do for [ii=1:STATS_records:1000] {
    # Jedes Objekt wird zweimal gezeichnet, einmal die Trajektorie und einmal ein Kreis am aktuellen Ort
    # Bei Erde und Mond wird der Kreis im korrekten Groessenverhaeltnis gemalt
    plot 'res/out2.dat' using 2:3 every ::1::ii w l lw 2 lc rgb "red" notitle, \
         'res/out2.dat' using 2:3 every ::ii::ii w p pt 6 ps 2 lc rgb "red" t "Sat", \
         'res/out2.dat' using 4:5 every ::1::ii w l lw 2 lc rgb "blue" notitle, \
         'res/out2.dat' using 4:5:(1.0) every ::ii::ii w circles lc rgb "blue" t "Earth", \
         'res/out2.dat' using 6:7 every ::1::ii w l lw 2 lc rgb "green" notitle,  \
         'res/out2.dat' using 6:7:(0.273) every ::ii::ii w circles lc rgb "green" t "Moon"
}
