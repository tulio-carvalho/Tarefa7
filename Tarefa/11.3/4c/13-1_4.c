//-y"+py'+qy+r=0 -> O método das diferenças finitas resolve equações desse tipo
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define N 19

double p( double x )
{	
	return -4/ x;
}

double q( double x )
{	
	return 2/ pow(x, 2);
}

double r( double x )
{	
	return -2*log(x)/ pow(x, 2);
}

double y( double x )
{	
	return 4/ x - 2/ pow(x, 2) + log(x) - 3./2.;
}

void **pivoteamento( double **M, int NL, int NC)
{
        double l, pivo, maior, aux;
        int i,j,k,m,n;
        
        for( j = 0; j < NL-1; j++)
        {	
		pivo = M[j][j];
		maior = pivo;

		for( k = j; k < NL; k++)
		{	if( fabs(maior) < fabs(M[k][j]))
			{	maior = M[k][j];
				m = k;
			}
		}

		if( maior != pivo)
		{	for( n = 0; n < NC; n++)
			{	aux = M[m][n];
				M[m][n] = M[j][n];
				M[j][n] = aux;
			}
			
		}

		for( i = j+1; i < NL; i++)
		{
			l = M[i][j]/M[j][j];       
			for( k = 0; k < NC; k++ )
				M[i][k] -= l*M[j][k]; 
		}
	}       
}

void subreversa( double **M, double *w, int dim, double *x)
{
	int i,j;
	double sum;

	w = malloc( dim* sizeof(double));
	for( i=dim; i>=0; i-- )
	{
		sum=0;
		for( j=i+1; j<dim; j++)
			sum += M[i][j]* w[j];

		w[i] = ( M[i][dim] - sum )/M[i][i];
	}
	//printf("x \t\ty\ty_exato\n");
	for( i = 1; i <= dim; i++ )
		printf("%lf\t%lf\t%lf\n",x[i],w[i-1],y(x[i]));
}

void main()
{	double	**M, x[N+2], h, w[N+2];
	int NL = N, NC = N+1, i, j;
	double xi = 1, xf = 2, yi = 1./2., yf = log(2);
 
	M = malloc( N* sizeof(double*));
	for( i = 0; i < NC; i++ )
		M[i] = malloc( NC* sizeof(double));
		
	h = 0.5;
	x[0] = xi;
	for( i = 1; i <= N; i++ )
	{	
		x[i] = xi + i*h;
	}
	x[N+1] = xf;
		
	//Criando a Matriz
	for( i = 0; i < N; i++ )
	{	
		if( i == 0 )
		{	M[i][i] = 2 + pow(h, 2)* q(x[i+1]);
			M[i][i+1] = - 1 + (h/2.)* p(x[i+1]);
		}
		else if( i == N-1 )
		{	M[i][i-1] = - 1 - (h/2.)* p(x[i+1]);
			M[i][i] = 2 + pow(h, 2)* q(x[i+1]);
		}
		else
		{	M[i][i-1] = - 1 - (h/2.)* p(x[i+1]);
			M[i][i] = 2 + pow(h, 2)* q(x[i+1]);
			M[i][i+1] = - 1 + (h/2.)* p(x[i+1]);
		}
	}
	M[0][N] = -pow(h, 2)* r(x[1]) + (1 +(h/2.)* p(x[1]))* y(x[0]);
	M[N-1][N] = -pow(h, 2)* r(x[N]) + (1 -(h/2.)* p(x[N]))* y(x[N+1]);
	
	for( i = 1; i < N-1; i++ )
	{	M[i][N] = -pow(h, 2)* r(x[i+1]);
		//printf("%lf\n", x[i]);
	}	

	pivoteamento(M,NL,NC);
	subreversa(M,w,NL,x);
}
