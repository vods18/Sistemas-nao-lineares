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
    double epsilon; //o epsilon para o criterios
    int max_iter; //numero maximo de iteracoes 
    //double** sl; //sistema lienar
    
}bag;

//void le_nome(int argc, char **argv, char* output);
void abre_arqs(FILE *arq, FILE*arq2);
void confere(FILE *arq, FILE *arq2);
double timestamp();
void clean_fgets(char *pos);
char* gera_incognitas(int max_eq, int w);
void cria_jacobiana(bag *b, char***jacobiana);
double* newton(bag *b, int cont_bag);
void trab1();
double norma_vetor(bag *b, double *x);
void analize_jacobiana_x(char*** jacobiana, double* x, char **names, int max_eq, double** values);
int split (const char *txt, char delim, char ***tokens);
void analize_function(bag *b, double *x, double *values, char **names, int cont_bag, int cont_aux);
double *eliminacaoGauss(bag *b, double** jacobiana_x, double *invert_x);
#endif // __UTILS_H__
