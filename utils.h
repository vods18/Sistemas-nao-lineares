#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdlib.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

typedef struct{
    int max_eq; //numero maximo de equacoes
    char** eq; //matriz de strings vinda com a entrada
    double* icognitas; //vetor com as icognitas
    char** jacobiana; //matriz jacobiana de derivadas
    double *x0; //iteracao inicial 
    double epsilon; //o epsilon para o criterio
    int max_iter; //numero maximo de iteracoes 
}bag;

double timestamp();
void clean_fgets(char *pos);
char** cria_jacobiana();
void newton (bag *b);
void trab1();
int split (const char *txt, char delim, char ***tokens);
#endif // __UTILS_H__
