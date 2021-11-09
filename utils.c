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

//calcula tempo de execucao em milisegundos
double timestamp(){
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}

//função para "limpar" string
void clean_fgets(char *pos) { 
  strtok(pos, "\n");
}

void gauss_pivo(){


}

void norma_vetor(){

    
}

char*** cria_jacobiana(bag *b){

  void *f, *f_dv;
  char ***jacob;
  double func;
  jacob = malloc(b->max_eq * sizeof(void*));
  for(int i=0; i<=b->max_eq; i++)
    jacob[i] = malloc(b->max_eq * sizeof(void));
   
  for(int i=0; i<b->max_eq; i++){
      clean_fgets(b->eq[i]);
      f = evaluator_create(b->eq[i]); //utilizamos as funções de cálculo de funções definidas pela biblioteca MATHEVAL
      assert(f);

      for(int j = 0; j < b->max_eq; j++){
        
        char var[100] = "x";  
        char *result;
        char num[100];
        int teste = j+1;
        sprintf(num, "%i", teste);
        strcat(var, num); //x1,x2,x3,...

        f_dv = evaluator_derivative (f, var); //também utilizamos essa biblioteca para calcular a derivada
        jacob[i][j] = evaluator_get_string(f_dv);
      }
  }
   
  return (jacob);
}

void newton (bag *b){

    for(int i=0; i<b->max_iter; i++){

        gauss_pivo(); 
        norma_vetor();

    }

}

void evaluator(bag *b){

    // double func;
    // double* linha;
    // for(int i=0; i<=b->max_eq; i++){
    //   clean_fgets(b->eq[i]);
    //   void *f = evaluator_create(b->eq[i]);
    //   assert(f);
    // }

    // for(int j=0; j<b->max_eq; j++){
    //   func = evaluator_evaluate_x(f, j);
    //   linha[j] = func;
    // }

}

int split (const char *txt, char delim, char ***tokens)
{
    int *tklen, *t, count = 1;
    char **arr, *p = (char *) txt;

    while (*p != '\0') if (*p++ == delim) count += 1;
    t = tklen = calloc (count, sizeof (int));
    for (p = (char *) txt; *p != '\0'; p++) *p == delim ? *t++ : (*t)++;
    *tokens = arr = malloc (count * sizeof (char *));
    t = tklen;
    p = *arr++ = calloc (*(t++) + 1, sizeof (char *));
    while (*txt != '\0')
    {
        if (*txt == delim)
        {
            p = *arr++ = calloc (*(t++) + 1, sizeof (char *));
            txt++;
        }
        else *p++ = *txt++;
    }
    free (tklen);
    return count;
}