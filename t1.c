#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include <inttypes.h>

int main (){

    //-- ler dados de sistemas.dat ----------------------------------------
    
    bag *b; //declaracao de ponteiro para a estrutura contendo variaveis de acordo com formato proposto

    //leitura das variáveis a partir de um arquivo 
    scanf("%i", &(b->max_eq)); //ler do arquivo dat maximo de equacoes possiveis

    for(int i=0; i<=b->max_eq; i++){

        char *equacao = malloc(sizeof(500)); //crio vetor auxiliar para ir recebendo por linha as funcoes dadas no dat
        fgets(equacao, 24, stdin);
        char ch;
        if(strlen(equacao) > 0){
        ch = equacao[0];
        }

        // analiza se foi feita a leitura de string inválida
        if(equacao == NULL || equacao == "" || equacao == " " || equacao == "\n" || equacao == "\0" || ch == 13 || ch==10){ 
        fgets(equacao, 24, stdin);
        }
        b->eq[i] = equacao;
    }

    for(int i=0; i<=b->max_eq; i++){
        scanf("%i", &(b->x0[i]));
    }

    
    scanf("%le", &(b->epsilon));
    scanf("%i", &(b->max_iter));

    printf("da");



    // evaluator();

    // --------------------------------------------------------------------

    // gerar matrix jacobiana com as equações fornecidas (libmatheval)
    // e->jacobiana = cria_jacobiana();
    // implementar Eliminação de Gauss com pivoteamento parcial. 
    // implementar função de cálculo da norma do vetor (do jeito que está no vídeo (12:26))

    // implementar  Método de Newton.
    // newton();
}