
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include<mpi.h>

#define NMAX 192   // Tamanho da malha
#define dMAX2 1/(NMAX*NMAX)
#define NMC 5000   // Numero de passos de Monte Carlo para cada T
#define dNMC 1/NMC
#define NP 192     // Numero de passos sobre o intervalo T 
#define SDOWN 0.1   // Quantidade de spin pra baixo 
#define ti 1.0	// Temperatura inicial
#define tf 5.0	// Temperatura final


//Variáveis para o processo de número aleatório
#define max 4294967295.0
#define dmax 1/max
unsigned s;	//Semente 
unsigned R;	//Raíz

//Variável para dar nomes aos arquivos do Gif 
int arq=1;


//Função que gera um número aleatório
double aleatorio(unsigned *R)  
{	double rn;

	*R=*R*s;
	rn=*R*dmax;

	return rn;
}

//Alocando a matrize
double ** Aloca_matriz(int linha, int coluna) 
{	int i;
        double **m;
        
        m = (double **)malloc(linha*sizeof(double*));
        for(i = 0; i < linha; i++)
                m[i] = (double *)malloc(coluna*sizeof(double));

        return m;
}

//Alocação o vetor
double * Aloca_vetor(int N)  
{	double *vet;
	
	vet = (double *)malloc(N*sizeof(double));
	
	return vet;
}
   


// Liberação de memória das matrizes e vetores utilizados
void Libera(double **Mat1, double **Mat2, double *v1, double *v2, double *v3, double *v4)
{	int i;

	for(i=0;i<(NMAX+2);i++) free(Mat1[i]);
	for(i=0;i<(NMAX+2);i++) free(Mat2[i]);
	free(Mat1);
	free(Mat2);
	free(v1);
	free(v2);
	free(v3);
	free(v4);
}

#include"ising_serial.h"

FILE *fp;

void main(int argc, char **argv)
{	int i, j;
	int rank, size, ierr, strip_size, root=0, to, from;
	double **spin, *spinD, **strip_S, *strip_D, *Vcima, *Vbaixo;
	double t, dt;
	MPI_Datatype strip;
	MPI_Status status;
	

	// Inicializando a região paralela
	ierr = MPI_Init(&argc,&argv);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	from=rank-1;
	if(from<0)
		from=size-1;
	
	to=rank+1;
	if(to>=size)
		to=0;
		
	//Raiz e Semente do número aleatório
	R = 893221891;
	s = (rank+1)*888121;
	
	Vcima = Aloca_vetor(NMAX);
	Vbaixo = Aloca_vetor(NMAX);
	
	if(rank==0)
	{	fp = fopen("saida.dat","w"); 
		strip_size = NMAX/size;
		
		spinD = Aloca_vetor(NMAX*(NMAX+2));  
		spin = (double **)malloc(sizeof(double*)*NMAX); 
		for(i=0;i<NMAX;i++)
		{	spin[i]=&(spinD[i*(NMAX+2)]);
		}
		// Inicialização
		inispin(spin); 		
	}
	
	MPI_Bcast(&strip_size,1,MPI_INT,0,MPI_COMM_WORLD);  
	MPI_Type_vector(strip_size,NMAX+2,NMAX+2,MPI_DOUBLE,&strip); 
	MPI_Type_commit(&strip);
	
	// Alocando vetor e matriz (stripdata e strip_S)
	strip_D = Aloca_vetor(strip_size*(NMAX+2));
	strip_S = (double **)malloc(sizeof(double*)*strip_size);
        for(i=0; i<strip_size; i++) {
                strip_S[i] = &(strip_D[i*(NMAX+2)]);
        }
        
        //Enviando com Scatter
	MPI_Scatter(spinD, 1, strip, &(strip_S[0][0]), 1, strip, 0, MPI_COMM_WORLD);

	
	dt=(tf-ti)/(NP-1);
	
	for(t=ti;t<=tf;t+=dt)
	{  	// Loop na quantidade de pontos gerado no gráfico
		MPI_Barrier(MPI_COMM_WORLD);
		MonteCarlo(strip_S,Vcima,Vbaixo,1/t,rank,size,from,to,strip_size,fp,&status);  
	}
	
	if(rank==0)
	{	MPI_Type_free(&strip);
		free(spin);
		free(strip_S);
		free(strip_D);
		fclose(fp);
		free(spinD);
		free(Vcima);
		free(Vbaixo);
	}
	
	MPI_Finalize(); 
}
