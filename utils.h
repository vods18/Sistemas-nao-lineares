#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdlib.h>
#include <inttypes.h>

typedef struct{
    int max_eq; //numero maximo de equacoes
    char** eq; //matriz de strings vinda com a entrada
    double* icognitas; //vetor com as icognitas
    char** jacobiana; //matriz jacobiana de derivadas
    double *x0; //iteracao inicial 
    double epsilon; //o epsilon1 para o criterio
    int max_iter; //numero maximo de iteracoes 
}bag;

double timestamp();
void clean_fgets(char *pos);
char** cria_jacobiana();
void newton (bag *b);
void trab1();

#endif // __UTILS_H__
