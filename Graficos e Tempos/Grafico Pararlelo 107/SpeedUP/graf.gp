set term pngcairo 
set output "MPISpeedUp107.png"

set xlabel "Core"
set ylabel "SpeedUp(Core)"
set xtics 1
set xrange [1:8]
set yrange [1:6]
set key right bottom

plot "MPI1Tempo107.dat" u 1:2 w l lw 2 lc "red" t "Temperatura", "MPI2Tempo107.dat" u 1:2 w l lw 2 lc "blue" t "Dominio"
