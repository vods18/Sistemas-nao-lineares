#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include <inttypes.h>
#include <assert.h>

#define MAX_BAG 10


int main (){

    //-- ler dados de sistemas.dat ----------------------------------------
    FILE *arq=fopen("sistemas.dat","r");

    if (!arq){
      perror ("Erro ao abrir arquivo sistemas.dat") ;
      exit (1) ; // encerra o programa com status 1
    }


    bag **bags; //vetor de bags, cada bag corresponde a um bloco da entrada
    bags = malloc (MAX_BAG * sizeof(bag*));
    
    for (int i=0; i<MAX_BAG; i++){
        bags[i] = malloc(sizeof(bag));
    }
    
    int cont_bag=1;

    int i = 0;
    while(!feof(arq)){
    
        bag *b = malloc(sizeof(bag)); //declaracao de ponteiro para a estrutura contendo variaveis de acordo com formato proposto

        // b->max_eq ----------------
        char max_eq[24];
        fgets(max_eq, 24, arq); 
        clean_fgets(max_eq);
        b->max_eq = atoi(max_eq);
        // --------------------------

        // b->eq ----------------------------------------------------------------------------------------------------------------------------
        b->eq = malloc (b->max_eq * sizeof(char*)); 
        for (int i=0; i<b->max_eq; i++){
            b->eq[i] = malloc(24 * sizeof(char));
        }

        for(int i=0; i<=b->max_eq -1; i++){
            char *equacao = malloc(sizeof(500)); //crio vetor auxiliar para ir recebendo por linha as funcoes dadas no dat
            fgets(equacao, 24, arq);
            char ch;

            if(strlen(equacao) > 0){
                ch = equacao[0];
            }


            // analiza se foi feita a leitura de string inválida
            if(equacao == NULL || equacao == "" || equacao == " " || equacao == "\n" || equacao == "\0" || ch == 13 || ch==10){ 
                fgets(equacao, 24, arq);
            }

            b->eq[i] = equacao;

        }
        // -------------------------------------------------------------------------------------------------------------------------------------


        // b->x0 -------------------------------------------------------------------------------------------------------------------------------
        b->x0 = malloc((b->max_eq -1) * sizeof(double)); 
        char* x0 = malloc(100 * sizeof(char));
        char ch;

        fgets(x0, 100, arq);

        if(x0 == NULL || x0 == "" || x0 == " " || x0 == "\n" || x0=="0" || x0 == "\0"){
            fgets(x0, 100, arq);
        }

        char **tokens;
        int count, i;
        const char *str = x0;

        count = split(str, ' ', &tokens);
        for (i = 0; i < count; i++){
            b->x0[i] = atof(tokens[i]);
        }
        
        // -------------------------------------------------------------------------------------------------------------------------------------

        // b->epsilon ----------------------------------------------------------------------------------------------------------------------------
        char ep[10];
        fgets(ep, 10, arq);
        b->epsilon = atof(ep);
        // -------------------------------------------------------------------------------------------------------------------------------------


        char max_iter[24];
        fgets(max_iter, 24, arq); //ler do arquivo dat maximo de equacoes possiveis
        clean_fgets(max_iter);
        b->max_iter = atoi(max_iter);

        bags[i] = b; 
        i++;

        printf("Bloco: %i\n", cont_bag);
        //--------teste entrada
        printf("Dimensão : %i\n", b->max_eq);
        for(int i=0; i<=b->max_eq -1; i++){
            printf("Equação %i: %s", i, b->eq[i]);
        }
        for(int i=0; i<=b->max_eq -1; i++){
            printf("x%i: %f\n", i, b->x0[i]);
        }

        printf("Epsilon: %f\n", b->epsilon);
        printf("Max_iter: %i\n", b->max_iter);
        
        char ext = fgetc(arq);

        free(b);
        cont_bag++;
        printf("\n");
    }

    // evaluator();

    // --------------------------------------------------------------------

    // gerar matrix jacobiana com as equações fornecidas (libmatheval)
    // e->jacobiana = cria_jacobiana();
    // implementar Eliminação de Gauss com pivoteamento parcial. 
    // implementar função de cálculo da norma do vetor (do jeito que está no vídeo (12:26))

    // implementar  Método de Newton.
    // newton();
}