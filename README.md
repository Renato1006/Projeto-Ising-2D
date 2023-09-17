<p>Projeto: Ising 2D</p>
<h1>Autor: Renato Maciel Félix</h1>

<p>Projeto final da matéria Computação de Alto Desempenho II (HPC2). Onde buscamos resolver um problema numérico e aplicar conceitos de paralelismo, com objetivo de obter maior performance no nosso código.</p>

<p>Problema: Comportamento de uma rede quadrada de spins, com a distribuição de 90% dos spins pra cima. É calculado a magnetização, energia e o calor específico, enquanto variamos a temperatura (ti=1, tf=5).</p>

<p>Máquina pessoal: Intel(R) Core(TM) i7-5500U CPU @ 2.40GHz</p>

<p>Após o uso do makefile, na pasta bin encontra-se o execultável. O resultado da simulação, com os valores das variaveis pela temperatura, também estará disponível lá (saida.dat).</p>

<p>Na pasta obj temos o nosso programa, onde em ising_serial.h temos uma bibioteca com as principais funções e em ising_serial.c temos o andamento da simulação.</p>

<p>Parte do projeto utiliza MPI para aplicação de paralelismo.</p>



<p>Foi realizado um teste, com uma rede de 100x100, com 1000 passos de Monte Carlo e uma divisão de 100 pontos, o programa levou 1m12.602s na máquina pessoal. Aumentando a rede para 256x256, com 5000 passos de Monte Carlo e com 250 pontos, o progrma levou 120m56.866s segundos na mesma máquina.</p>

<p>O projeto encontra-se concluido, podemos ver o resultado observando o pdf junto ao diretório.</p>