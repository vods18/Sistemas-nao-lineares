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

void cria_jacobiana(bag *b, char***jacobiana){ //OK

  void *f, *f_dv;
  double func;
   
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
      f_dv = evaluator_derivative(f, var); //também utilizamos essa biblioteca para calcular a derivada
      assert(evaluator_get_string(f_dv));
      jacobiana[i][j] = evaluator_get_string(f_dv);
    }
    evaluator_destroy(f);
  }
}


void analize_function(bag *b, double *x, double *values, char **names, int cont_bag, int cont_aux){ // TODO: ERRO TA AQUI

  double val=0, maior = 0;
  printf("Contador interno: %i\n", cont_aux);
  // substitui nas equações originais os valores de x0
  void *f;
  for(int i = 0; i< b->max_eq; i++){
    // if(i == 1){
    //   // void  *g = malloc(24  * sizeof(void));
    //   // printf("id g: %x\n", g);
    //   // clean_fgets(b->eq[i]);
    //   // printf("-------->%s\n", b->eq[i]);
    //   // g = evaluator_create(b->eq[i]); /*SEGFAULT*/ 
    //   // printf("UFA\n");
    //   // assert(g);
    //   // val = evaluator_evaluate(g, b->max_eq, names , x); 
    //   // values[i] = val;
    //   // evaluator_destroy(g);
    // } else{

    clean_fgets(b->eq[i]);
    // printf("id f: %x\n", f);
    // printf("%s\n", b->eq[i]);
    f = evaluator_create(b->eq[i]);
    assert(f);
    val = evaluator_evaluate(f, b->max_eq, names , x); 
    values[i] = val;
    evaluator_destroy(f);
   //}  
  }
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

  // printf("\n--------- names ------------------------------=--------------\n");
  // for(int h = 0; h < max_eq; h++){
  //   printf("%s\n", names[h]);
  // }
  // printf("\n-------------------------------------------------------------\n");

  // printf("\n--------- Jacobiana ---------------------------------------\n");
  // for(int i =0; i< max_eq; i++){
  //   for(int j =0; j < max_eq; j++){
  //     printf("%s      ", jacobiana[i][j]);
  //   }
  //   printf("\n");
  // }
  // printf("-------------------------------------------------------------\n"); 

  for(int i =0; i< max_eq; i++){
    for(int j =0; j < max_eq; j++){
      clean_fgets(jacobiana[i][j]);
      void *f;
      double val;
      // assert(jacobiana[i][j]);
      // printf("i= %d , j = %d\n", i , j);
      assert(jacobiana[i][j]);
      f = evaluator_create(jacobiana[i][j]); //utilizamos as funções de cálculo de funções definidas pela biblioteca MATHEVAL
      val = evaluator_evaluate(f, max_eq, names , x); 
      values[i][j] = val; 
      evaluator_destroy(f);
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

    free(matrix);
    free(vetorB);
    return x;
}

double* newton (bag *b, FILE* arq2, int cont_bag){
    int cont_aux=0;
    void *f;
    double *x = b->x0; // valor calculado na iteração anterior x1,x2,x3,... (x0)
    double *delta = malloc((b->max_eq -1) * sizeof(double)); // valor calculado na iteração atual para x1,x2,x3,...
    double *x_novo = malloc((b->max_eq -1) * sizeof(double)); // x + delta
    double *values = malloc((b->max_eq -1) * sizeof(double));
    double *invert_x = malloc((b->max_eq -1) * sizeof(double));
    double **jacobiana_x = malloc((b->max_eq -1) * sizeof(double*));    
    for(int s=0; s< b->max_eq; s++){
      jacobiana_x[s] = malloc((b->max_eq -1) * sizeof(double));
    }
   
    char **incognitas = malloc(MAX_LIN * sizeof(char*));  
    for(int j=0; j< b->max_eq; j++){
      incognitas[j] = malloc(MAX_LIN * sizeof(char));
    }

    char ***jacobina  = (char ***) malloc(sizeof(char**) * b->max_eq); //FIX 1
    for(int i = 0; i < b->max_eq; i++){
      jacobina[i] = (char **) malloc(sizeof(char*) * b->max_eq);  // FIX 2
      for(int j = 0; j < b->max_eq; j++){
        jacobina[i][j] = (char *) malloc(sizeof(char) * 100);
      }
    }
    
    b->tderivadas = timestamp();
    cria_jacobiana(b, jacobina);
    b->tderivadas = timestamp() - b->tderivadas;
    
    for(int i=0; i<b->max_iter; i++){

      printf("#\n");
      int inter=1;
      for(int s=0; s< b->max_eq; s++){
        printf("x%d = %f\n", inter , x[s]);
        inter++;
      }
      printf("#\n");

      // incognitas = [x1, x2, x3, ..]
      for(int w=0; w<b->max_eq; w++){
        char var[MAX_LIN] = "x";  
        char num[MAX_LIN];
        int teste = w+1;
        sprintf(num, "%i", teste);
        strcat(var, num);
        for(int z=0; z<3; z++){
          incognitas[w][z]=var[z];
        }
      }
      // analize_function(b,x, values, incognitas, cont_bag, cont_aux); //f(x)

      // ---------------------------------------------
      double val=0, maior = 0;
      // printf("Contador interno: %i\n", cont_aux);
      // substitui nas equações originais os valores de x0
      if(cont_bag <3){
        for(int pt = 0; pt< b->max_eq; pt++){
          clean_fgets(b->eq[pt]);
          // printf("id f: %x\n", f);
          // printf("%s\n", b->eq[pt]);
          f = evaluator_create(b->eq[pt]);
          assert(f);
          val = evaluator_evaluate(f, b->max_eq, incognitas , x); 
          values[pt] = val;
          evaluator_destroy(f);
        }
      } else {
        void *g = (void*) malloc(sizeof(void*));
        for(int pq = 0; pq< b->max_eq; pq++){
          clean_fgets(b->eq[pq]);
          // printf("id f: %x\n", f);
          // printf("%s\n", b->eq[pq]);
          g = evaluator_create(b->eq[pq]);
          assert(g);
          val = evaluator_evaluate(g, b->max_eq, incognitas , x); 
          values[pq] = val;
          evaluator_destroy(g);
        }
      }
      // ---------------------------------------------

      cont_aux++;
      // printf("------------>SAÍ EIN\n");

      if(norma_vetor(b, values) < b->epsilon){
        puts("A");
        return x;
      }

      // jacobiana(x) 
      double tjacobina = timestamp();
      analize_jacobiana_x(jacobina, x, incognitas, b->max_eq, jacobiana_x);   
      b->tjacobiana = b->tjacobiana + (timestamp() - tjacobina);
      
      // -f(x)
      for(int m = 0; m< b->max_eq; m++){
        invert_x[m] = ((-1) * values[m]);
      }

      // jacobiana(x) * incognitas = - f(x) => SL

      double tsl = timestamp();
      delta = eliminacaoGauss(b, jacobiana_x, invert_x);
      b->tsl = b->tsl + (timestamp() - tsl);
        
      for(int a = 0; a<b->max_eq; a++){
        x_novo[a] = delta[a] + x[a];
      }

      if(norma_vetor(b, delta)< b->epsilon){
        puts("B");
        return x_novo;
      }

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
    free(delta);
    free(x_novo);
    free(values);
    free(invert_x);
    for(int s=0; s< b->max_eq; s++){
      free(jacobiana_x[s]);
    }
    free(jacobiana_x);
    /*for(int i = 0; i < b->max_eq; i++){
      free(jacobina[i]);
      for(int j = 0; j < b->max_eq; j++){
        free(jacobina[i][j]);
      }
    }*/
    for(int j=0; j< b->max_eq; j++){
      free(incognitas[j]);
    }
    
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

