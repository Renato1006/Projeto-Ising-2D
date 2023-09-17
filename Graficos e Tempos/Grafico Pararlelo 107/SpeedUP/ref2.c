#include<stdio.h>
int main()
{	double ref, x, y;
	int a, b;
	
	do
	{
		printf("\nValor ref: ");
		scanf("%lf",&ref);

		do
		{
			printf("\nValor otimizado: ");
			scanf("%lf",&x);
 
			y = ref/x; 
			

			printf("\nSpeedUP: %lf\n", y);
			
		}while(x!=0);
		
	}while(ref!=0);
}
