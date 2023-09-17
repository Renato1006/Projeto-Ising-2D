
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
int arq=0;


//Função que gera um número aleatório
double aleatorio(unsigned *R)  
{	double rn;

	*R=*R*s;
	rn=*R*dmax;

	return rn;
}

//Alocando a matrize
double ** Aloca_matriz() 
{	int i;
        double **m;
        
        m = (double **)malloc((NMAX+2)*sizeof(double*));
        for(i = 0; i < (NMAX+2); i++)
                m[i] = (double *)malloc((NMAX+2)*sizeof(double));

        return m;
}

//Alocação o vetor
double * Aloca_vetor()  
{	double *vet;
	
	vet = (double *)malloc(NP*sizeof(double));
	
	return vet;
}

// Liberação de memória das matrizes e vetores utilizados
void Libera(double **Mat1, double **Mat2, double *v1, double *v2, double *v3, double *v4, double *v5)
{	int i;

	for(i=0;i<(NMAX+2);i++) free(Mat1[i]);
	for(i=0;i<(NMAX+2);i++) free(Mat2[i]);
	free(Mat1);
	free(Mat2);
	free(v1);
	free(v2);
	free(v3);
	free(v4);
	free(v5);
}

#include"ising_serial.h"

FILE *fp;

void main(int argc, char **argv)
{	int i, j, k, rank, size, ierr, count, root=0;
	double **spin, **Fut, *mag, *energ, *Temp, *C, *Tt;
	double dt, T, temp, x, m_atual, e_atual;
	
	//Raiz e Semente do número aleatório
	R = 893221891;
	s = 888121;
	
	// Inicializando a região paralela
	ierr = MPI_Init(&argc,&argv);
        ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
        ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
        s = 888121*(rank+1);
	
	fp = fopen("saida.dat","a");
	
	count = NP/size;

	//Alocando memória para as matrizes e vetores do problema
	spin=Aloca_matriz();
	Fut=Aloca_matriz();		
	mag=Aloca_vetor(count);
	energ=Aloca_vetor(count);
	C=Aloca_vetor(count);
	Temp=Aloca_vetor(NP);
	Tt=Aloca_vetor(count);

	inispin(spin); 

	temp=ti;
	dt=(tf-ti)/(NP-1);
	
	for(i=0;i<NP;i++)
	{	Temp[i]=1/temp;
		temp+=dt;
	}
	
	//Enviando com Scatter
	ierr = MPI_Scatter(Temp,count,MPI_DOUBLE,Tt,count,MPI_DOUBLE,root,MPI_COMM_WORLD);
	
	arq=1+rank*count;
	
	for(i=0;i<count;i++) 
	{	MonteCarlo(spin,Fut,mag,energ,Tt[i],i,C,&arq);  
		fprintf(fp,"%lf\t%lf\t%lf\t%lf\n", 1/Tt[i], mag[i], energ[i], C[i]); 
	}
	
	// Finalizando a região paralela
	ierr = MPI_Finalize();
	
	Libera(spin,Fut,mag,energ,C,Temp,Tt);
	
	fclose(fp);
}
