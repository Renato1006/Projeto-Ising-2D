set term pngcairo 
set output "EfiGrafFlags.png"

set xlabel "OX"
set ylabel "Eficiência (%)"
set xtics 1
set key right bottom

plot "EfiTempoFlagOX.dat" u 1:2 w l lw 2 lc "red" t "OX", "EfiTempoFlagOXAvx.dat" u 1:2 w l lc "blue" lw 2 t "Avx", "EfiTempoFlagOXNative.dat" u 1:2 w l lc "green" lw 2 t "Native"
