#!/bin/bash

#mpicc -o test1.x ising_serial.c -O3 -fexpensive-optimizations -m64 -foptimize-register-move -funroll-loops -ffast-math -mtune=native -march=native -lm


export FI_SOCKETS_IFACE=eth0
export FI_PROVIDER=tcp

mkdir TesteUltimo

for i in 1 2 3 4 5 6 7 8 
do
	echo "Rodando com $i cores."
	time mpirun -np $i ./ising.x
	mkdir Saida$i
	mv saida.dat Saida$i
	mv Saida$i TesteUltimo
	echo "" 
done

cp nohup.out TesteUltimo
