set term pngcairo 
set output "MPItempo107.png"

set xlabel "Core"
set ylabel "Tempo(Segundos)"
set xtics 1

plot "MPI1Tempo107.dat" u 1:2 w l lw 2 lc "red" t "Temperatura", "MPI2Tempo107.dat" u 1:2 w l lw 2 lc "blue" t "Dominio"
