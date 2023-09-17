set term pngcairo 
set output "EfiMPItempo107.png"

set xlabel "Core"
set ylabel "Eficiência (%)"
set xtics 1
set key right bottom
set yrange [0:100]

plot "EfiMPI1Tempo107.dat" u 1:2 w l lw 2 lc "red" t "Temperatura", "EfiMPI2Tempo107.dat" u 1:2 w l lw 2 lc "blue" t "Dominio"
