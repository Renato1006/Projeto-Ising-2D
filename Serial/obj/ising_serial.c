
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>

#define NMAX 192  // Tamanho da malha
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

void main()
{	int i, j;
	double **spin, **Fut, *mag, *energ, *Temp, *C;
	double dt, T, temp, x, m_atual, e_atual;
	
	//Raiz e Semente do número aleatório
	R = 893221891;
	s = 888121;
	
	fp = fopen("saida.dat","w"); // Abrindo o arquivo de saída

	//Alocando memória para as matrizes e vetores do problema
	spin=Aloca_matriz();
	Fut=Aloca_matriz();		
	mag=Aloca_vetor();
	energ=Aloca_vetor();
	C=Aloca_vetor();
	Temp=Aloca_vetor();

	// Inicialização
	inispin(spin); 

	T=ti;
	dt=(tf-ti)/(NP-1); // Infinitésimo de temperatura a ser adicionado a cada passo
	
	for(i=0;i<NP;i++) 
	{	temp=1/T; //Passando 1/T como parametro de MonteCarlo()
		MonteCarlo(spin,Fut,mag,energ,temp,i,C,&arq);  //Calculos do programa
		fprintf(fp,"%lf\t%lf\t%lf\t%lf\n", T, mag[i], energ[i], C[i]); //Salvando as variaveis
		T+=dt; //Atualizando a temperatura com dT
	}
	
	Libera(spin,Fut,mag,energ,C,Temp);
	
	fclose(fp);
}
