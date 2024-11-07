#include <omp.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <climits> // INT_MAX

using namespace std;

class Celula {
    public:
        int i, j;
    
        Celula () {
            
        }
        
        bool operator==(Celula cel) {
            return (this->i == cel.i && this->j == cel.j);
        }
        
        bool operator!=(Celula cel) {
            return (this->i != cel.i || this->j != cel.j);
        }
};

typedef struct fila {
    Celula cel;
    fila *prox;  
} fila;


fila* verificaVizinhos(Celula *celEtapAtual, int **Grid, int m, int n);

int main(int argc, char** argv) {
    if(argc != 3) {
        printf("O programa foi executado com parâmetros incorretos.\n");
        printf("Uso: ./rotpar <arquivo_de_entrada> <arquivo_de_saída>\n");
        exit(1);
    }
    
    ifstream entrada(argv[1]);
    ofstream saida;
    
    // Células de origem e de destino
    Celula origem, destino;

    int m, n, i, j;

    // Leitura do tamanho da matriz
    entrada >> m >> n;
    
    int **Grid = new int*[m];

    // Leitura dos índices da origem
    entrada >> origem.i >> origem.j;
    
    // Leitura dos índices do destino
    entrada >> destino.i >> destino.j;
    
    int q;
    
    // Leitura da quantidade de obstáculos
    entrada >> q;
    
    int obstI[q], obstJ[q], qtdLin[q], qtdCol[q];
    
    for (int c = 0; c < q; c++) {
        entrada >> obstI[c] >> obstJ[c] >> qtdLin[c] >> qtdCol[c];
    }

    // Variáveis utilizadas para medir o tempo de execução do programa
    double tini, tfin, tseq;
    
    tini = omp_get_wtime();
    
    // Criação das linhas da matriz
    for (i = 0; i < m; i++) {
        Grid[i] = new int[n];
    }

    #pragma omp parallel
    {
        // Inicializa todas as células com valor 'infinito'
        #pragma omp for private(j)
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++) {
                Grid[i][j] = INT_MAX;
            }
        
        // Inicialização dos q obstáculos
        #pragma omp for private(i, j)
        for (int c = 0; c < q; c++) {

            // Inicializa os obstáculos com -1
            for (i = 0; i < qtdLin[c]; i++) {
                for (j = 0; j < qtdCol[c]; j++) {
                    Grid[obstI[c] + i][obstJ[c] + j] = -1;
                }
            }
        }
    }

    
    /**** Busca pela célula de destino ****/
    
    // Inicializa a distância da origem
    Grid[origem.i][origem.j] = 0;
   
    fila *ini = new fila;
    ini->prox = NULL;
    fila *fim = ini;
    ini->cel = origem;
    fila *p;
    
    int k = 0;  // Etapa atual: cada etapa corresponde a distância já percorrida da origem até o destino
    int qtdCelEtapa; // Quantidade de células da etapa atual presentes na fila

    const int SIZE = 300000;
    fila *resultado[SIZE] = {NULL};
    Celula celEtapAtual[SIZE];
    
    bool achou = false;
    int etapa_atual;

    // Expansão das células    
    while (ini != NULL && !achou) {
        p = ini;
        qtdCelEtapa = 0;
        etapa_atual = Grid[p->cel.i][p->cel.j];
        
        // Se p != NULL, então ele aponta para alguma célula, porém é preciso verificar
        // se ele aponta para uma célula pertencente a esta etapa (k),
        // e se a "quantidade de células desta etapa" (qtdCelEtapa) é menor que o tamanho máximo do vetor.
        // Ou seja, se tem espaço sobrando no vetor celEtapAtual. Pois caso não tenha, ela(s) será(rão) adicionada(s) na(s) próxima(s) iteração(ões) do
        // while mais externo: while (ini != NULL && !achou)
        while (p != NULL && k == etapa_atual && qtdCelEtapa < SIZE) {
            celEtapAtual[qtdCelEtapa].i = p->cel.i;
            celEtapAtual[qtdCelEtapa].j = p->cel.j;
            qtdCelEtapa++;
            p = p->prox;
        }
        
        /* Verifico se posso incrementar o contador de etapas, para indicar que todas as células da etapa atual já estão a
         * caminho de serem expandidas pelo "for", logo abaixo.
        */
        // Se p continha apenas células da etapa atual
        if (p == NULL) {
            k++;
          // Senão, se apesar de p conter células de etapas posteriores, todas as células da etapa atual já foram adicionadas em celEtapAtual,
          // pois atualmente ele aponta para uma célula da etapa seguinte
        } else if (k != Grid[p->cel.i][p->cel.j]) {
            k++;
        }
        
        #pragma omp parallel for
        for (int i = 0; i < qtdCelEtapa; i++) {
            resultado[i] = verificaVizinhos(&celEtapAtual[i], Grid, m, n);  // retorna subfilas contendo os vizinhos de cada célula analisada
        }
            
        int i = 0;
        // Junta todas as subfilas em uma só
        while (i < qtdCelEtapa) {
            fim->prox = resultado[i];
            if (fim->prox != NULL) {
                fim = fim->prox;
                if (fim->prox != NULL) {
                    fim = fim->prox;
                    if (fim->prox != NULL) {
                        fim = fim->prox;
                    }
                }
            }
            i++;
        }
            
        // Remove os primeiros elementos da fila que já tiveram seus vizinhos verificados (expandidos)
        if (ini != NULL) {
            while (qtdCelEtapa-- > 0) {
                p = ini;
                ini = ini->prox;
                
                if (p->cel == destino) {
                    achou = true;
                }
                delete p;
            }
        }
    };
    
    // Caso tenha encontrado
    if (achou) {
        // Libera memória ocupada pela fila
        for ( ; ini != NULL; ) {
            p = ini;
            ini = ini->prox;
            
            delete p;
        }

        
        /**** Backtracking ****/
        
        ini = new fila;
        ini->prox = NULL;
        fim = ini;
        ini->cel = destino;
        
        // Vetor que guardará o caminho de volta à origem (caso exista), que consiste na quantidade de passos até o destino + 1
        fila caminho[Grid[destino.i][destino.j] + 1];
        
        caminho[Grid[destino.i][destino.j]].cel = destino;  // Guarda o destino na última posição de caminho[].cel
        i = Grid[destino.i][destino.j];
        p = ini;
        
        int caminhoI, caminhoJ;
        
        while (i > 0) {
                caminhoI = caminho[i].cel.i;
                caminhoJ = caminho[i].cel.j;

                // Se não for o extremo superior da matriz, ou seja, se não for a primeira linha da matriz E seu valor for Infinito
                if (caminhoI-1 >= 0 && Grid[caminhoI-1][caminhoJ] == Grid[caminhoI][caminhoJ] - 1) {
                    caminho[i-1].cel.i = caminho[i].cel.i-1;
                    caminho[i-1].cel.j = caminho[i].cel.j;
                }
                
                // Se não for o extremo direito da matriz, ou seja, se não for a última coluna da matriz E seu valor for Infinito
                else if (caminhoJ+1 < n && Grid[caminhoI][caminhoJ+1] == Grid[caminhoI][caminhoJ] - 1) {
                    caminho[i-1].cel.i = caminho[i].cel.i;
                    caminho[i-1].cel.j = caminho[i].cel.j+1;
                }
                
                // Se não for o extremo inferior da matriz, ou seja, se não for a última linha da matriz E seu valor for Infinito
                else if (caminhoI+1 < m && Grid[caminhoI+1][caminhoJ] == Grid[caminhoI][caminhoJ] - 1) {
                    caminho[i-1].cel.i = caminho[i].cel.i+1;
                    caminho[i-1].cel.j = caminho[i].cel.j;
                }
                
                // Se não for o extremo esquerdo da matriz, ou seja, se não for a primeira coluna da matriz E seu valor for Infinito
                else if (caminhoJ-1 >= 0 && Grid[caminhoI][caminhoJ-1] == Grid[caminhoI][caminhoJ] - 1) {
                    caminho[i-1].cel.i = caminho[i].cel.i;
                    caminho[i-1].cel.j = caminho[i].cel.j-1;
                }
                
                i--;
        };

        tfin = omp_get_wtime();
        tseq = tfin - tini;
        printf("Tempo total: %f\n\n", tseq);
        
        // Impressão do backtracking no arquivo de saída
        saida.open(argv[2], std::ofstream::trunc);
        saida << Grid[destino.i][destino.j] << endl;
        
        for (i = 0; i <= Grid[destino.i][destino.j]; i++) {
            saida << caminho[i].cel.i << " " << caminho[i].cel.j << endl;
        }
        saida.close();
    }
    
    if (!achou) {
        tfin = omp_get_wtime();
        tseq = tfin - tini;
        printf("Não existe um caminho até o destino.\nTempo total: %f\n\n", tseq);
    }
    
    // Libera memória ocupada pela matriz
    for(i = 0; i < m; i++) {
        delete [] Grid[i];
    }
    delete [] Grid;
    
    return 0;
}

fila* verificaVizinhos(Celula *celEtapAtual, int **Grid, int m, int n) {
    fila *ini = NULL;
    fila *fim = NULL;
    
    // Se existir uma linha acima, ou seja, se houver uma célula E seu valor for Infinito
    if (celEtapAtual->i-1 >= 0 && Grid[celEtapAtual->i-1][celEtapAtual->j] == INT_MAX) {
        Grid[celEtapAtual->i-1][celEtapAtual->j] = Grid[celEtapAtual->i][celEtapAtual->j] + 1;
        ini = new fila;
        fim = ini;

        fim->cel.i = celEtapAtual->i-1;
        fim->cel.j = celEtapAtual->j;
    }

    // Se existir uma coluna à esquerda, ou seja, se houver uma célula E seu valor for Infinito
    if (celEtapAtual->j-1 >= 0 && Grid[celEtapAtual->i][celEtapAtual->j-1] == INT_MAX) {
        Grid[celEtapAtual->i][celEtapAtual->j-1] = Grid[celEtapAtual->i][celEtapAtual->j] + 1;
        if (ini == NULL) {
            ini = new fila;
            fim = ini;
        } else {
            fim->prox = new fila;
            fim = fim->prox; 
        }
        
        fim->cel.i = celEtapAtual->i;
        fim->cel.j = celEtapAtual->j-1;
    }

    // Se existir uma coluna à direita, ou seja, se houver uma célula E seu valor for Infinito
    if (celEtapAtual->j+1 < n && Grid[celEtapAtual->i][celEtapAtual->j+1] == INT_MAX) {
        Grid[celEtapAtual->i][celEtapAtual->j+1] = Grid[celEtapAtual->i][celEtapAtual->j] + 1;
        if (ini == NULL) {
            ini = new fila;
            fim = ini;
        } else {
            fim->prox = new fila;
            fim = fim->prox; 
        }
        
        fim->cel.i = celEtapAtual->i;
        fim->cel.j = celEtapAtual->j+1;
    }

    // Se existir uma linha abaixo, ou seja, se houver uma célula E seu valor for Infinito
    if (celEtapAtual->i+1 < m && Grid[celEtapAtual->i+1][celEtapAtual->j] == INT_MAX) {
        Grid[celEtapAtual->i+1][celEtapAtual->j] = Grid[celEtapAtual->i][celEtapAtual->j] + 1;
        if (ini == NULL) {
            ini = new fila;
            fim = ini;
        } else {
            fim->prox = new fila;
            fim = fim->prox; 
        }
        
        fim->cel.i = celEtapAtual->i+1;
        fim->cel.j = celEtapAtual->j;
    }
    
    // Caso a célula atual tenha algum vizinho
    if (fim != NULL)
        fim->prox = NULL;

    return ini;
}