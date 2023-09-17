set term pngcairo 
set output "GrafFlags.png"

set xlabel "OX"
set ylabel "Tempo(Segundos)"
set xtics 1

plot "TempoFlagOX.dat" u 1:2 w l lw 2 lc "red" t "OX", "TempoFlagOXAvx.dat" u 1:2 w l lc "blue" lw 2 t "Avx", "TempoFlagOXNative.dat" u 1:2 w l lc "green" lw 2 t "Native"
