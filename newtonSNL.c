#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include <inttypes.h>
#include <assert.h>

#define MAX_NOME 100 //maximo de caracteres para um nome de arquivo de saida

int main (int argc, char **argv){

    //ler dados de sistemas.dat -------------------------------------------------------------------------------------------------------------
    FILE *arq, *arq2;
    
    //abre_arqs(arq,arq2);
    arq = fopen("sistemas.dat","r");

    char* output = malloc(MAX_NOME * sizeof(char));
    output=malloc(MAX_NOME * sizeof(char)); // reservo espaço para um nome de ate 100 letras
    output = le_nome(argc, argv);

    if (output == NULL){
        puts("SOU NULL");
		arq2 = stdout; //caso nao tenha sido passado um nome, pegue da saida padrao
    } else{
        puts(output);
        arq2 = fopen(output, "w"); //Crio arquivo
    }

    int cont_bag = 0;
    while(!feof(arq)){
        bag *b = malloc(sizeof(bag)); //declaracao de ponteiro para a estrutura contendo variaveis de acordo com formato proposto

        b->ttotal = 0;
        b->ttotal = timestamp();


        // b->max_eq ------------------------------------------------------------------------------------------------------------------------
        char max_eq[24];
        fgets(max_eq, 24, arq); 
        clean_fgets(max_eq);
        b->max_eq = atoi(max_eq);
        // ----------------------------------------------------------------------------------------------------------------------------------

        // b->eq ----------------------------------------------------------------------------------------------------------------------------
        b->eq = malloc (b->max_eq * sizeof(char*)); 
        for (int i=0; i<b->max_eq; i++){
            b->eq[i] = malloc(24 * sizeof(char));
        }

        for(int i=0; i<=b->max_eq -1; i++){
            b->eq[i] = malloc(sizeof(500)); //crio vetor auxiliar para ir recebendo por linha as funcoes dadas no dat
            fgets(b->eq[i], 24, arq);
            char ch;

            if(strlen(b->eq[i]) > 0){
                ch = b->eq[i][0];
            }


            // analiza se foi feita a leitura de string inválida
            if(b->eq[i] == NULL || b->eq[i] == "" || b->eq[i] == " " || b->eq[i] == "\n" || b->eq[i] == "\0" || ch == 13 || ch==10){ 
                fgets(b->eq[i], 24, arq);
            }

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

        // b->epsilon --------------------------------------------------------------------------------------------------------------------------
        char ep[10];
        fgets(ep, 10, arq);
        b->epsilon = atof(ep);
        // -------------------------------------------------------------------------------------------------------------------------------------

        // b->max_iter -------------------------------------------------------------------------------------------------------------------------
        char max_iter[24];
        fgets(max_iter, 24, arq); //ler do arquivo dat maximo de equacoes possiveis
        clean_fgets(max_iter);
        b->max_iter = atoi(max_iter);
        // -------------------------------------------------------------------------------------------------------------------------------------


        // Prints ------------------------------------------------------------------------------------------------------------------------------
        fprintf(arq2,"\n---------Início do bloco ------------------------------------\n");
        fprintf(arq2,"Dimensao:  %d\n", b->max_eq);
        for(int h = 0; h < b->max_eq; h++){
            fprintf(arq2,"Equacao %i : %s\n", h + 1, b->eq[h]);
        }
        // -------------------------------------------------------------------------------------------------------------------------------------

        // Método de Newton. -------------------------------------------------------------------------------------------------------------------
        if(cont_bag < 3) {
            newton(b, arq2, cont_bag);
        }
        
        b->ttotal = timestamp() - b->ttotal;

        // -------------------------------------------------------------------------------------------------------------------------------------
        // Print tempos ------------------------------------------------------------------------------------------------------------------------
        fprintf(arq2, "###########\n");
        fprintf(arq2, "# Tempo Total: %f\n", b->ttotal);
        fprintf(arq2, "# Tempo Derivadas: %f\n", b->tderivadas);
        fprintf(arq2, "# Tempo Jacobiana: %f\n", b->tjacobiana);
        fprintf(arq2, "# Tempo SL: %f\n", b->tsl);   
        fprintf(arq2, "###########\n");     
        fprintf(arq2, "-------------------------------------------------------------\n");
        // -------------------------------------------------------------------------------------------------------------------------------------

        char ext = fgetc(arq);
        for(int c=0; c<b->max_eq; c++)
            free(b->eq[c]);
        free(b);
        free(x0);
        cont_bag++; 
    }
    // -----------------------------------------------------------------------------------------------------------------------------------------

    fclose(arq);
    fclose(arq2);
}