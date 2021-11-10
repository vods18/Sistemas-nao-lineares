#include "utils.h"
#include <string.h>
#include <math.h>
#include <matheval.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include <sys/time.h>
#include <stdlib.h>

// double* gauss_pivo(bag *b,double** jacobiana_x, char **incognitas, double *invert_x){
int eliminacaoGauss (bag *b, double** jacobiana_x, char **incognitas, double *invert_x){
double **matrix = jacobiana_x;
double *vetorB = invert_x;

int i,j,k,d,e;
unsigned int n;
n = b->max_eq;    matrix[k][i] = 0.0;

    for( int j=i+1; j < n; ++j ){
        matrix[k][j] -= matrix[i][j] * m;
    }

    vetorB[k] -= vetorB[i] * m;

    }
}


// Calculates x from x[n-1] to x[0]
for (int i = n - 1; i >= 0; --i) {
    real_t s = 0;
    for (int j = i + 1; j < n; ++j){
        s = s + matrix[i][j]*x[j];
    }
    x[i] = (vetorB[i]-s)/matrix[i][i];
}

double temp = timestamp() - *tTotal;
printf("===> Eliminação Gauss: %.17g ms\n",temp );
printf("     -->X : ");
prnVetor(x, sl->n);
printf("\n");

return 0;
}