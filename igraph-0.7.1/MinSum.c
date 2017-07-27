/* 
Algoritmo Min-Sum:
dif_i = vecindario del vertice i
I = funcion indicadora -> 1 si el argumento es vdd y 0 sino.
S es un conjunto yde deciclicaje de G si X*_i(S) = 1 para todo vertice i.
Xt_i son variables binarias.  t_i(S) = min(t : Xt_i(S) = 1)

De hecho, en el procedimiento de extracción de las hojas, se retira un vértice 
en el primer paso después del momento en que todos los vecinos, excepto uno, 
han sido retirados, convirtiéndolo en una hoja.

partition function -> se computa a partir del intercambio de mensajes entre los nodos
vecinos. Esos mensajes son interpretados como las leyes de probabilidad marginal
de las variables locales donde algunas interacciones deben ser removidas.
Los mensajes entrantes de un nodo se tratan como independientes y la energia libre se
puede calcular como una suma de contribuciones locales dependiendo de la solucion del mensaje
n_ij(t_i,t_j) es el mensaje de i a j 

h_i(t_i) es el costo factible minimo de remover nodos que hacen que el 2-core del
grafo resultante sea vacio.
El nodo i pertenece al conjunto de deciclicaje optimo si cumple ecc 9

Min-Sum computa los valores aproximados de h_i mediante la resoluciones de un sistema
de ecuaciones  -> quedan las ecc 10a y 10b. Para esto se debe resolver el sistema de ecuaciones
dado por las ecc 11a - 11g
Darse un T arbitrario tq 0 < t_i <= T. phi_i(t_i) = I[t_i = 0](...)+...

La solucion del sistema se obtiene mediante iteraciones. Es suficiente tomar un T <= N
*/
#include <igraph.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void print_vector(igraph_vector_t *v, FILE *f) {
  long int i;
  for (i=0; i<igraph_vector_size(v); i++) {
    fprintf(f, " %li", (long int) VECTOR(*v)[i]);
  }
  fprintf(f, "\n");
}

/* Se encarga de calcular el tamaño de la componente conexa mas grande */
int max_component(igraph_t *g){
	igraph_vector_ptr_t complist;

	/* Se inicializa vector que contendra todas las componentes conexas del grafo */
	igraph_vector_ptr_init(&complist, 0);
	igraph_decompose(g,&complist,IGRAPH_WEAK,-1,2); // se realiza la descomposicion del grafo en sus componentes

	/* Se verifica cual es la componente mas grande */
	int max = 0;
	for(int i=1; i<igraph_vector_ptr_size(&complist); i++){
		if(igraph_vcount(VECTOR(complist)[max]) < igraph_vcount(VECTOR(complist)[i])){ // Hay una mas grande que la componente actual
			max = i;
		}
	}
	return (int)(igraph_vcount(VECTOR(complist)[max])); // componente conexa mas grande
}


int main(){
	FILE *F;
	char filename[32];
	igraph_t graph;
	igraph_vector_t degrees, nodes, CIvalues;
	double remove = 0.1; // multiplicador de porcentaje
	double rem_nodes = 0.0; // cantidad de nodos removidos
	int total_nodes; // total de nodos del grafo original
	int T = 35; // limite de la iteracion
	clock_t start, end;
	double time_used;

	G = fopen("MinSum_times.csv","w"); // archivo que guardara los tiempo de ejecucion por iteracion
	H = fopen("MinSum_iter.csv", "w"); // archivo que guardara los R-index (componente mas grande) por iteracion

	/* Se lee el archivo que contiene las conexiones de los nodos */
	F = fopen("red3.edges","r");
	igraph_read_graph_edgelist(&graph,F,0,0); // crea el grafo a partir del archivo con las conexiones
	fclose(F);

	total_nodes = igraph_vcount(&graph);

	/* Proceso de escritura del nuevo grafo tras cierto porcentaje de eliminacion de nodos */
		rem_nodes += 1.0; // aumento cantidad de nodos removidos
		if(rem_nodes == ceil(total_nodes*remove)){
			sprintf(filename,"grafo%d_MinSum.edges",(total_nodes - (int)rem_nodes)); // nombre del archivo donde estaran los resultados
		    F = fopen(filename,"w");
		    igraph_write_graph_edgelist(&graph,F); // escritura
		    fclose(F);
		    remove += 0.1; // aumento porcentaje de eliminacion
		}



	
	fprintf(stderr, "%i %i\n", igraph_vcount(&graph), igraph_ecount(&graph));

	printf("VACIO\n");

	return 0;
}