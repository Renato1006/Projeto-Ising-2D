//Calculo da soma dos spin dos vizinhos   
double vizinho(double **S, int linha, int coluna) 
{	double v;

	v = S[linha+1][coluna]+S[linha-1][coluna]+S[linha][coluna+1]+S[linha][coluna-1];
	
	return v;
}

// Inicialização da matriz
void inispin(double **S)  
{	int i, j;

	for(i=1;i<=NMAX;i++)
        {       for(j=1;j<=NMAX;j++)
                {      //Definindo a posição inicial dos spins.
			if(aleatorio(&R)>=SDOWN) 
                      		S[i][j]=1;
               	else
                       	S[i][j]=-1;  
                }
	}

	// Inicializando os elementos da matriz estendida
	for(i=1;i<=NMAX;i++)
	{	S[i][0] = S[i][NMAX];
		S[i][NMAX+1] = S[i][1];
	}
	
	for(j=1;j<=NMAX;j++)
	{	S[0][j] = S[NMAX][j];
		S[NMAX+1][j] = S[1][j];
	}
	
	S[0][0] = 0;
	S[0][NMAX+1] = 0;
	S[NMAX+1][0] = 0;
	S[NMAX+1][NMAX+1] = 0;
}

//Inicialização da magnetização e energia para cada interação
void inicia_int(double **S, double *m, double *e) 
{	int i,j;
	
	*m=0;
	*e=0;

	for(i=1;i<=NMAX;i++)
	{	for(j=1;j<=NMAX;j++)
		{	// Calculo da magnetização
			*m+=S[i][j];  
			// Cálculo da energia
			*e=(-1)*((S[i][j]*vizinho(S,i,j))*0.5); 
		}
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
void Calor_Especifico(double *Calor, int l, double Ea, double Ea2, double T)
{	double  x;

	//Variação da energia
	x=(Ea2*dNMC)-((Ea*dNMC)*(Ea*dNMC));

	Calor[l]=(x*T*T)*dMAX2;  
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
void MonteCarlo(double **S, double **Fut, double *m, double *e, double t, int l, double *calor, int *c)
{	int i, j, k;
	double cont=NMC*NMAX*NMAX, dcont=1/cont;
	double ma, ea, m_sum=0, e_sum=0, e_sum2=0, soma, dE;

	for(i=0;i<NMC;i++)
	{	//Calculo da energia e magnetização para cada nova iteração
		inicia_int(S,&ma,&ea);  

		//Preparando a matriz Futuro
		for(k=0;k<(NMAX+2);k++)
			for(j=0;j<(NMAX+2);j++)
				Fut[k][j]=S[k][j];
				
 		// Primeira linha		
		for(k=1;k<=NMAX;k++)
		{	soma=vizinho(Fut,1,k);
                       dE=Fut[1][k]*soma*2;

                       Metropolis(Fut,1,k,&ma,&ea,dE,t);
		}
		
		// Atualização das linhas extendidas
		for(j=1;j<=NMAX;j++)
		{	Fut[NMAX+1][j]=Fut[1][j];
			Fut[0][j]=Fut[NMAX][j];
		}
		
		// Meio da matriz
		for(j=2;j<=NMAX-1;j++)
		{	for(k=1;k<=NMAX;k++)
			{	soma=vizinho(Fut,j,k);
				dE=Fut[j][k]*soma*2;

				Metropolis(Fut,j,k,&ma,&ea,dE,t);
			}
			
			// Atualização das colunas extendidas com os extremos de cada linha
			Fut[j][0]=Fut[j][NMAX];
			Fut[j][NMAX+1]=Fut[j][1];
		}

		// Ultima linha
		for(k=1;k<=NMAX;k++)
                {	soma=vizinho(Fut,NMAX,k);
                       dE=Fut[NMAX][k]*soma*2;

                       Metropolis(Fut,NMAX,k,&ma,&ea,dE,t);
                }
		
		// Atualização das linhas extendidas
		for(j=1;j<=NMAX;j++)
		{	Fut[NMAX+1][j]=Fut[1][j];
			Fut[0][j]=Fut[NMAX][j];
		}

		//Variáveis para o cálculo da média e do calor específico
             	m_sum+=ma; 
		e_sum+=ea;
		e_sum2+=(ea)*(ea); 
	}
	
	//Atualizando a matriz S com as atualizações da mtariz Futuro
	for(i=0;i<(NMAX+2);i++)
		for(j=0;j<(NMAX+2);j++)
			S[i][j]=Fut[i][j];
	
	//Salvando o spin do microestado
	//salvaSpin(S,*&c); 

	
	//Encontrando a magnetização e energia média
	m[l]= m_sum*dcont; 
	e[l]= e_sum*dcont;

	Calor_Especifico(calor, l, e_sum, e_sum2, t); // Calculando o calor específico para a T atual
}
