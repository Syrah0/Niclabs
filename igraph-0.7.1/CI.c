/* 

Algoritmo CIl: Debe tener l inicial -> consideraremos l = 4 (por comparacion realizada en paper de COREHD)
1.- Calcular el CI value de los nodos -> calcular CI nodo i a partir de sus vecinos de distancia l.
2.- Construir max-heap
3.- Remover la raiz (la idea es remover del grafo el nodo con mayor CIl).
4.- Heapify -> Recalcular los valores de CI de los nodos -> Heapify y volver a 3.

*/
#include <igraph.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void print_vector(igraph_vector_t *v, FILE *f) {
  long int i;
  for (i=0; i<igraph_vector_size(v); i++) {
    fprintf(f, " %li", (long int) VECTOR(*v)[i]);
  }
  fprintf(f, "\n");
}

/* Se encarga de seleccionar los nodos que se deben eliminar para conformar el 2-core */
igraph_vector_t nodesToDistance(igraph_t *g, int l, int node){
	igraph_vector_t nodes, nodesExc, neigh;
	//igraph_vector_ptr_t neighborhood;

	/* inicializa el vector que contendra los nodos a distancia l */
	igraph_vector_init(&nodes, 0);
	igraph_vector_push_back(&nodes, node);

	/* inicializa el vector que contendra los nodos ya visitados */
	igraph_vector_init(&nodesExc, 0);
	igraph_vector_push_back(&nodesExc, node); // ingresa nodo a calcular el CI

	/* revisar resto de nodos */
	for(int i = 1; i <= l; i++){
		igraph_vector_t nodesaux;
		/* inicializa el vector con nodos del vecindario de node */
		igraph_vector_init(&nodesaux, 0);

		/* calcula el vecindario del nodo node a distancia i */
		for(int j = 0; j < igraph_vector_size(&nodes); j++){
			igraph_vector_init(&neigh, 0);
			igraph_neighbors(g, &neigh, igraph_vector_e(&nodes, j), IGRAPH_ALL);
			for(int k = 0; k < igraph_vector_size(&nodesExc); k++){
				for(int m = 0; m < igraph_vector_size(&neigh); m++){
					if(igraph_vector_e(&neigh,m) == igraph_vector_e(&nodesExc, k)){
						igraph_vector_remove(&neigh, m);
					}
				}
			}
			for(int k = 0; k < igraph_vector_size(&neigh); k++){
				igraph_vector_push_back(&nodesExc, igraph_vector_e(&neigh, k));
				igraph_vector_push_back(&nodesaux, igraph_vector_e(&neigh, k));
			}
			igraph_vector_destroy(&neigh);
		}

		igraph_vector_destroy(&nodes); // destruye vector para realizar copia
		igraph_vector_copy(&nodes, &nodesaux);

		igraph_vector_destroy(&nodesaux);
		igraph_vector_destroy(&neigh);
	}

	igraph_vector_destroy(&nodesExc);
	return nodes;
}

/* Se encarga de calcular la componente conexa mas grande */
igraph_t* component(igraph_t *g){
	igraph_vector_ptr_t complist;

	/* Se inicializa vector que contendra todas las componentes conexas del grafo */
	igraph_vector_ptr_init(&complist, 1);
	igraph_decompose(g,&complist,IGRAPH_WEAK,-1,2); // se realiza la descomposicion del grafo en sus componentes

	/* Se verifica cual es la componente mas grande */
	int max = 0;
	for(int i=1; i<igraph_vector_ptr_size(&complist); i++){
		if(igraph_vector_size(VECTOR(complist)[max]) < igraph_vector_size(VECTOR(complist)[i])){ // Hay una mas grande que la componente actual
			max = i;
		}
	}
	return VECTOR(complist)[max]; // componente conexa mas grande
}


int main(){
	FILE *F;
	igraph_t graph;
	igraph_vector_t degrees, nodes, CIvalues;
	int l = 2;

	/* CIl = (Ki - 1) * SUM_{j vecinos de i a distancia l} (Kj - 1) */
	/* Se lee el archivo que contiene las conexiones de los nodos */
	F = fopen("red2.edges","r");
	igraph_read_graph_edgelist(&graph,F,0,0); // crea el grafo a partir del archivo con las conexiones
	fclose(F);

	/* Calcular los grados de los nodos del grafo */
	igraph_vector_init(&degrees, 0);
	igraph_degree(&graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS); 

	/* Calcular el CI de cada nodo */
	igraph_vector_init(&CIvalues, 0);
	for(int i = 0; i < igraph_vector_size(&degrees); i++){
		nodes = nodesToDistance(&graph, l, i);
		int CI = (igraph_vector_e(&degrees, i) - 1);
		int sum = 0;
		for(int j = 0; j < igraph_vector_size(&nodes); j++){
			sum += (igraph_vector_e(&degrees, igraph_vector_e(&nodes, j)) - 1);
		}
		CI *= sum;
		igraph_vector_push_back(&CIvalues, CI);
	}
	print_vector(&CIvalues, stdout);

	printf("VACIO\n");

	/* realizar tree-breaking o desmantelamiento */
	return 0;
}