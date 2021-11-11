#ifndef __UTILS_H__
#define __UTILS_H__

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
#include <unistd.h>


typedef struct{
    int max_eq; //numero maximo de equacoes
    char** eq; //matriz de strings vinda com a entrada
    double* icognitas; //vetor com as icognitas
    char** jacobiana; //matriz jacobiana de derivadas
    double *x0; //iteracao inicial 
    double epsilon; //o epsilon para o criterios
    int max_iter; //numero maximo de iteracoes 
    double ttotal; //tempo total
    double tderivadas; //tempo das derivadas
    double tjacobiana; //tempo da jacobiana
    double tsl; //tempo da SL
    
}bag;

char *le_nome(int argc, char **argv);
void confere(FILE *arq, FILE *arq2);
double timestamp();
void clean_fgets(char *pos);
char* gera_incognitas(int max_eq, int w);
void cria_jacobiana(bag *b, char***jacobiana);
double* newton(bag *b, FILE*arq2, int cont_bag);
void trab1();
double norma_vetor(bag *b, double *x);
void analize_jacobiana_x(char*** jacobiana, double* x, char **names, int max_eq, double** values);
int split (const char *txt, char delim, char ***tokens);
void analize_function(bag *b, double *x, double *values, char **names, int cont_bag, int cont_aux);
double *eliminacaoGauss(bag *b, double** jacobiana_x, double *invert_x);
#endif // __UTILS_H__
