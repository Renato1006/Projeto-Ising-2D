
//Calculo da soma dos spin dos vizinhos   
double vizinho(double **S, int linha, int coluna) 
{	double v;

	v = S[linha+1][coluna]+S[linha-1][coluna]+S[linha][coluna+1]+S[linha][coluna-1];
	
	return v;
}

// Inicialização da matriz
void inispin(double **S)  
{	int i, j;

	for(i=0;i<NMAX;i++)
        {       for(j=1;j<=NMAX;j++)
                {      //Definindo a posição inicial dos spins.
			if(aleatorio(&R)>=SDOWN) 
                      		S[i][j]=1;
               	else
                       	S[i][j]=-1;  
                }
	}

	// Inicializando os elementos da matriz estendida
	for(i=0;i<NMAX;i++)
	{	S[i][0] = S[i][NMAX];
		S[i][NMAX+1] = S[i][1];
	}
	/*
	for(j=1;j<=NMAX;j++)
	{	S[0][j] = S[NMAX][j];
		S[NMAX+1][j] = S[1][j];
	}
	*/
	
	/*S[0][0] = 0;
	S[0][NMAX+1] = 0;
	S[NMAX+1][0] = 0;
	S[NMAX+1][NMAX+1] = 0;
	*/
}

//Inicialização da magnetização e energia para cada interação
void inicia_int(double **S, double *m, double *e, int size, double *vc, double *vb) 
{	int i,j,k;
	double soma=0;
	
	*m=0;
	*e=0;

	for(k=1;k<=NMAX;k++)
	{	*m+=S[0][k];
		soma=S[0][k-1]+S[0][k+1]+S[1][k]+vc[k-1];
		*e=(-1)*(S[0][k]*soma*0.5);
	}
	
	for(i=1;i<size-1;i++)
	{	for(j=1;j<=NMAX;j++)
		{	// Calculo da magnetização
			*m+=S[i][j];  
			// Cálculo da energia
			*e=(-1)*((S[i][j]*vizinho(S,i,j))*0.5); 
		}
	}
	
	// Ultima linha  (USA Vb)
	for(k=1;k<=NMAX;k++)
	{	*m+=S[size-1][k];
		soma=S[size-1][k-1]+S[size-1][k+1]+S[size-2][k]+vb[k-1];
		*e=(-1)*(S[size-1][k]*soma*0.5);
	}
}

// Gerar arquivos com a malha de spins para CADA microestado (ou cada temperatura)
void salvaSpin(double **S, int *k)
{	int i,j;
	FILE *fp2;
	char arq[30];

	sprintf(arq,"arq%d.dat", *k);  // gerando arq1.dat, arq2.dat ......
	fp2 = fopen(arq,"w");

	for(i=1;i<=NMAX;i++) // Escrevendo a matriz de spin no arquivo correspondente
        {       for(j=1;j<=NMAX;j++)
        	        fprintf(fp2,"%.0lf\t", S[i][j]);
                fprintf(fp2,"\n");
        }
	*k=*k+1;

	fclose(fp2);	
}

// Calculo do calor especifico 
double Calor_EspecificoP(double Ea, double Ea2, double T)
{	double  x, C;

	//Variação da energia
	x=(Ea2*dNMC)-((Ea*dNMC)*(Ea*dNMC));

	C = (x*T*T)*dMAX2;
	
	return C;  
}

// Algoritmo de Metropolis para atualização das variáveis de interesse
void Metropolis(double **M, int linha, int coluna, double *m, double *e, double dE, double t)
{	double prob;
	
	if(dE <= 0)
	{       M[linha][coluna] = -M[linha][coluna];
                *m += M[linha][coluna]*2;
                *e -= dE;
	}
        else
        {      //Cálculo do fator de boltzmann 
        	prob = exp((-1)*(dE*t)); 
        
                if(prob>aleatorio(&R))
                {       M[linha][coluna] = -M[linha][coluna];
                        *m += M[linha][coluna]*2;
                        *e += dE;
                }
        }
}

//Função de Monte Carlo
void MonteCarlo(double **S, double *Vc, double *Vb, double t, int rank, int size, int from, int to, int Ssize, FILE *fp, MPI_Status *stat)
{	int i, j, k;
	double cont=NMC*NMAX*NMAX, dcont=1/cont;
	double ma, ea, m_sum=0, e_sum=0, e_sum2=0, soma, dE;
	double mtot=0, etot=0, m, e, c=0, *v1, *v2, **Fut, cv=0;
	
	//Alocando auxiliares
	v1 = Aloca_vetor(NMAX);
	v2 = Aloca_vetor(NMAX);
	Fut = Aloca_matriz(Ssize,NMAX+2);

	for(i=0;i<NMC;i++)
	{	//Preparando a matriz Futuro
		for(k=0;k<Ssize;k++)
			for(j=0;j<(NMAX+2);j++)
				Fut[k][j]=S[k][j];
		
		// Vetores vx recebem os elementos das colunas estendidas
 		for(k=1;k<=NMAX;k++)
 		{	v1[k-1]=S[0][k];
			v2[k-1]=S[Ssize-1][k];
		}

		// Enviando as tiras
		MPI_Send(v1,NMAX,MPI_DOUBLE,to,9,MPI_COMM_WORLD);
		MPI_Send(v2,NMAX,MPI_DOUBLE,to,99,MPI_COMM_WORLD);
		
		//Recebendo as tiras
		MPI_Recv(Vc,NMAX,MPI_DOUBLE,from,9,MPI_COMM_WORLD,stat);
		MPI_Recv(Vb,NMAX,MPI_DOUBLE,from,99,MPI_COMM_WORLD,stat);
	
		//Calculo da energia e magnetização para cada nova iteração
		inicia_int(S,&ma,&ea,Ssize,Vc,Vb);  

				
 		//Primeira linha		
		for(k=1;k<=NMAX-1;k++)
		{	soma = Fut[0][k-1] + Fut[0][k+1] + Fut[1][k] + Vc[k-1];
			dE=Fut[0][k]*soma*2;

			Metropolis(Fut,0,k,&ma,&ea,dE,t);
		}
		
		//Extremos
		Fut[0][NMAX+1] = Fut[0][1];
		soma = Fut[0][NMAX-1] + Fut[0][NMAX+1] + Fut[1][NMAX] + Vc[NMAX-1];
		dE = Fut[0][NMAX]*soma*2;
		Metropolis(Fut,0,NMAX,&ma,&ea,dE,t);
		
		//Meio da matriz
		for(j=1; j<Ssize-1; j++)
		{	for(k=1; k<=NMAX-1; k++)
			{	dE=Fut[j][k]*vizinho(Fut,j,k)*2;

				Metropolis(Fut,j,k,&ma,&ea,dE,t);
			}
			
			//Atualização das colunas extendidas com os extremos de cada linha
			soma = vizinho(Fut,j,NMAX);
			Fut[j][NMAX+1] = Fut[j][1];
			dE = Fut[j][NMAX]*soma*2;
			Metropolis(Fut,j,NMAX,&ma,&ea,dE,t);
		}

		//Ultima linha
		for(k=1; k<=NMAX; k++)
                {	soma = Fut[Ssize-1][k-1] + Fut[Ssize-1][k+1] + Fut[Ssize-2][k] + Vb[k-1];
                       dE = Fut[Ssize-1][k-1]*soma*2;

                       Metropolis(Fut,Ssize-1,k,&ma,&ea,dE,t);
                }
                
                //Extremos
                Fut[Ssize-1][NMAX+1] = Fut[Ssize-1][1];
                soma = Fut[Ssize-1][NMAX-1] + Fut[Ssize-1][NMAX+1] + Fut[Ssize-2][NMAX] + Vb[NMAX-1];
                
                dE = Fut[Ssize-1][NMAX]*soma*2;
                Metropolis(Fut,Ssize-1,NMAX,&ma,&ea,dE,t);
		
		//Variáveis para o cálculo da média e do calor específico
             	m_sum += ma; 
		e_sum += ea;
		e_sum2 += (ea)*(ea); 
	}
	
	c = Calor_EspecificoP(e_sum, e_sum2, t);
	
	MPI_Barrier(MPI_COMM_WORLD);

	// Reduce nas variaveis calculadas para o root
	MPI_Reduce(&m_sum,&mtot,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&e_sum,&etot,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&c,&cv,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	//Atualizando a matriz S com as atualizações da mtariz Futuro
	for(i=0;i<Ssize;i++)
		for(j=0;j<(NMAX+2);j++)
			S[i][j]=Fut[i][j];
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	//Salvando o spin do microestado
	//salvaSpin(S,*&c); 
	
	//Encontrando a magnetização e energia média
	if(rank==0)
	{	m = mtot*dcont;
		e = etot*dcont;
		
		//Gravando as variáveis
		fprintf(fp,"%lf\t%lf\t%lf\t%lf\n", 1/t, m, e, cv); 
	}

	free(v1);
	free(v2);
	for(i=0;i<Ssize;i++)
		free(Fut[i]);
	free(Fut);
}

