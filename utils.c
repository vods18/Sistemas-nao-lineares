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

#define MAX_NOME 20 //maximo de caracteres para um nome de arquivo de saida
#define MAX_LIN 100 //maximo de caracteres permitidas na linha de uma funcao passada no .dat

/*
void le_nome(int argc, char **argv, char* output){
  
  int option;
  
  while ((option = getopt (argc, argv, "o:")) != -1) {		

    switch(option) { 

      case 'o':

        output = optarg; //pego o nome da imagem para escrita
    }
  }
  
}*/

void abre_arqs(FILE *arq, FILE*arq2){
  
  char* output = malloc(MAX_NOME * sizeof(char));
  output=malloc(MAX_NOME * sizeof(char)); // reservo espaço para um nome de ate 30 letras
	output = "sample.out";
  //output = le_nome(argc, argv);

  arq = fopen("sistemas.dat","r");

  /*if (output == NULL)
		arq = stdout; //caso nao tenha sido passado um nome, pegue da saida padrao
	else*/
		arq2 = fopen(output, "w"); //Crio arquivo

}

void confere(FILE *arq, FILE *arq2){
    if (!arq){
      perror ("Erro ao abrir arquivo de entrada") ;
      exit (1) ; // encerra o programa com status 1
    }

    if (!arq2){
      perror ("Erro ao criar arquivo de saida") ;
      exit (1) ; // encerra o programa com status 1
    }
}

//calcula tempo de execucao em milisegundos
double timestamp(){ //OK
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}

//função para "limpar" string
void clean_fgets(char *pos) { //OK
  strtok(pos, "\n");
}

char*** cria_jacobiana(bag *b){ //OK

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
        
        char var[MAX_LIN] = "x";  
        char num[MAX_LIN];
        int teste = j+1;
        sprintf(num, "%i", teste);
        strcat(var, num); //x1,x2,x3,...

        f_dv = evaluator_derivative (f, var); //também utilizamos essa biblioteca para calcular a derivada
        jacob[i][j] = evaluator_get_string(f_dv);
      }
  }
   
  return (jacob);
}


double* analize_function(bag *b, double *x){ //OK

  char **eq = b->eq;
  void *f, *f_dv;
  double *values = malloc(b->max_eq * sizeof(double));
  double val=0, maior = 0;
  int i;
  
  char *name = malloc(2 * sizeof(char)); 
  char **names = malloc(MAX_LIN * sizeof(char*));  
  for(int j=0; j< b->max_eq; j++){
    names[j] = malloc(MAX_LIN * sizeof(char));
  }

  for(int w=0; w<b->max_eq; w++){
    char *name = malloc(2 * sizeof(char)); 
    char var[MAX_LIN] = "x";  
    char num[MAX_LIN];
    int teste = w+1;
    sprintf(num, "%i", teste);
    strcat(var, num);
    name = var;
    for(int z=0; z<2; z++){
      names[w][z]=name[z];
    }
  }

  // substitui nas equações originais os valores de x0
  for(i =0; i< b->max_eq; i++){
    clean_fgets(b->eq[i]);
    f = evaluator_create(b->eq[i]); //utilizamos as funções de cálculo de funções definidas pela biblioteca MATHEVAL
    assert(f);
    val = evaluator_evaluate(f, b->max_eq, names , x); 
    values[i] = val; 
  } 

  return values;
}

double norma_vetor(bag *b, double *x){ //OK
  double maior = 0;
  for(int j=0; j<b->max_eq; j++){
    if(fabs(x[j])> maior){
      maior = fabs(x[j]);
    }
  }
  return maior;
}

void analize_jacobiana_x(char*** jacobiana, double* x, char **names, int max_eq, double** values){ //OK

  for(int i =0; i< max_eq; i++){
    for(int j =0; j < max_eq; j++){
      clean_fgets(jacobiana[i][j]);
      void *f;
      double val;
      f = evaluator_create(jacobiana[i][j]); //utilizamos as funções de cálculo de funções definidas pela biblioteca MATHEVAL
      assert(f);
      val = evaluator_evaluate(f, max_eq, names , x); 
      values[i][j] = val; 
    }
  }
}

double *eliminacaoGauss(bag *b, double** jacobiana_x, double *invert_x){

    double **matrix = malloc((b->max_eq - 1) * sizeof(double*));
    for(int j=0; j< b->max_eq; j++){
      matrix[j] = malloc((b->max_eq - 1) * sizeof(double));
    }

    double *vetorB = malloc((b->max_eq - 1) * sizeof(double));
    double *x = malloc((b->max_eq - 1) * sizeof(double));

    int i,j,k,d,e;
    unsigned int n;
    n = b->max_eq;

    for(int i=0; i < n; ++i ) {
      for(int k=0; k < n; ++k ){
        matrix[i,k] = jacobiana_x[i,k];
      }
      vetorB[i] = invert_x[i];
    }

    // printf("\n--------- matrix(x) ---------------------------------------\n");
    // for(int h = 0; h < b->max_eq; h++){
    //   for(int r = 0; r < b->max_eq; r++){
    //     printf("%le    ", matrix[h][r]);
    //   }
    //   printf("\n");
    // }
    // printf("-------------------------------------------------------------\n"); 
    // printf("\n--------- vetorB --------------------------------------------\n");
    // for(int h = 0; h < b->max_eq; h++){
    //   printf("%le ", vetorB[h]);
    // }
    // printf("\n-------------------------------------------------------------\n");


    // pivoteamento parcial---------------------------
    for(int i=0; i < n; ++i ) {
      for(int k=i+1; k < n; ++k ) {
        unsigned int iPivo = 0;
        for(int d=i+1; d<n ; d++){
          if (abs(matrix[d][i]) > abs(matrix[i][i])){
            iPivo = d;
          }
        }
        if (i < iPivo){
          double aux;
          for(int r = 0; r < n; r++){
              aux = matrix[i][r];
              matrix[i][r] = matrix[iPivo][r];
              matrix[iPivo][r] = aux;
          }

          // Troca os elementos do vetor: b
          aux = vetorB[i];
          vetorB[i] = vetorB[iPivo];
          vetorB[iPivo] = aux;
        }
      }
    }
    //-------------------------------------------------
 
    for(int i=0; i < n; ++i ) {
      for(int k=i+1; k < n; ++k ) {

        double m = matrix[k][i] / matrix[i][i];
        matrix[k][i] = 0.0;

        for( int j=i+1; j < n; ++j ){
          matrix[k][j] =  matrix[k][j] - matrix[i][j] * m;
        }

        vetorB[k] -= vetorB[i] * m;

        // printf("\n--------- matrix(x) ---------------------------------------\n");
        // for(int h = 0; h < b->max_eq; h++){
        //   for(int r = 0; r < b->max_eq; r++){
        //     printf("%le    ", matrix[h][r]);
        //   }
        //   printf("\n");
        // }
        // printf("-------------------------------------------------------------\n"); 
      }
    }


    // Calculates x from x[n-1] to x[0]
    for (int i = n - 1; i >= 0; --i) {
      double s = 0;
      for (int j = i + 1; j < n; ++j){
          s = s + matrix[i][j]*x[j];
      }
      x[i] = (vetorB[i]-s)/matrix[i][i];
    }

    return x;
}

double* newton (bag *b, char*** jacobiana){

    double *x = b->x0; // valor calculado na iteração anterior x1,x2,x3,... (x0)
    double *delta = malloc((b->max_eq -1) * sizeof(double)); // valor calculado na iteração atual para x1,x2,x3,...
    double *x_novo = malloc((b->max_eq -1) * sizeof(double)); // x + delta
    double *values = malloc((b->max_eq -1) * sizeof(double));
    double *invert_x = malloc((b->max_eq -1) * sizeof(double));
    double **jacobiana_x = malloc((b->max_eq -1) * sizeof(double*));    
    for(int s=0; s< b->max_eq; s++){
      jacobiana_x[s] = malloc((b->max_eq -1) * sizeof(double));
    }
   

    char *name = malloc(2 * sizeof(char)); 
    char **incognitas = malloc(MAX_LIN * sizeof(char*));  
    for(int j=0; j< b->max_eq; j++){
      incognitas[j] = malloc(MAX_LIN * sizeof(char));
    }
    
    for(int i=0; i<b->max_iter; i++){

      printf("#\n");
      int inter=1;
      for(int s=0; s< b->max_eq; s++){
        printf("x%d = %f\n", inter , x[s]);
        inter++;
      }
      printf("#\n");

      values = analize_function(b,x); //f(x)

      // printf("\n--------- X -----------------------------------------------\n");
      // for(int h = 0; h < b->max_eq; h++){
      //   printf("%le ", values[h]);
      // }
      // printf("\n-------------------------------------------------------------\n");

      if(norma_vetor(b, values) < b->epsilon){
        return x;
      }

      // incognitas = [x1, x2, x3, ..]
      for(int w=0; w<b->max_eq; w++){
        char *name = malloc(2 * sizeof(char)); 
        char var[MAX_LIN] = "x";  
        char num[MAX_LIN];
        int teste = w+1;
        sprintf(num, "%i", teste);
        strcat(var, num);
        name = var;
        for(int z=0; z<2; z++){
          incognitas[w][z]=name[z];
        }
      }

      // jacobiana(x) 
      analize_jacobiana_x(jacobiana, x, incognitas, b->max_eq, jacobiana_x);   
      
      // printf("\n--------- Jacobiana(x) --------------------------------------\n");
      // for(int h = 0; h < b->max_eq; h++){
      //   for(int r = 0; r < b->max_eq; r++){
      //     printf("%le    ", jacobiana_x[h][r]);
      //   }
      //   printf("\n");
      // }
      // printf("-------------------------------------------------------------\n"); 
      
      // -f(x)
      for(int m = 0; m< b->max_eq; m++){
        invert_x[m] = ((-1) * values[m]);
      }
      // printf("\n--------- invert_x ------------------------------------------\n");
      // for(int h = 0; h < b->max_eq; h++){
      //   printf("%le ", invert_x[h]);
      // }
      // printf("\n-------------------------------------------------------------\n");

      // jacobiana(x) * incognitas = - f(x) => SL

      delta = eliminacaoGauss(b, jacobiana_x, invert_x);
      // printf("\n--------- Delta ---------------------------------------------\n");
      // for(int h = 0; h < b->max_eq; h++){
      //   printf("%le ", delta[h]);
      // }
      // printf("\n-------------------------------------------------------------\n");
        
      for(int a = 0; a<b->max_eq; a++){
        x_novo[a] = delta[a] + x[a];
      }

      // printf("\n--------- x/_novo -------------------------------------------\n");
      // for(int h = 0; h < b->max_eq; h++){
      //   printf("%le ", x_novo[h]);
      // }
      // printf("\n-------------------------------------------------------------\n");

      if(norma_vetor(b, delta)< b->epsilon){
        return x_novo;
      }

      // delta  = [-1.625   ,   -1.375] 
      // x_novo = [-0.625   ,   -3.625]

      for(int f=0; f<b->max_eq; f++)
        x[f] = x_novo[f];

    }

    printf("#\n");
    int inter=1;
    for(int s=0; s< b->max_eq; s++){
      printf("x%d = %f\n", inter , x[s]);
      inter++;
    }
    printf("#\n");
    return x;
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

