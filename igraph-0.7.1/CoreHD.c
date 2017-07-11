/* 

1-core : todos los nodos
2-core : remover todos los nodos con un arco o menos. Remover arcos asociados a ese nodo.

Algoritmo CoreHD:
1.- Obtener el 2-core -> obtener grados de cada nodo del 2-core
2.- Obtener nodo con mayor grado
3.- Remover dicho nodo -> actualizar 2-core y grados de sus nodos -> verificar nvo 2-core: Vacío -> tree-breaking (dismantling)
							   															   No vacío -> Volver a paso 2	

*/

#include <igraph.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(){
	FILE *F;
	igraph_t graphaux, graph;
	igraph_vector_t result, resultaux, remaux;
	igraph_vs_t rem;

	F = fopen("red3.edges","r");
	igraph_read_graph_edgelist(&graph,F,0,0); // crea el grafo a partir del archivo con las conexiones
	fclose(F);

	igraph_vector_init(&resultaux, 0);
	
	/* calculo de los grados de cada nodo del grafo */
	igraph_degree(&graph, &resultaux, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS); 

	int removed = 0; // servira para llevar la cuenta de la cantidad de nodos a remover
	//printf("Degree is  %10i, vertex %2i.\n", (int)igraph_vector_e(&resultaux,393),393);

	/* Revisa la cantidad de nodos de remover para formar el 2-core */
	for(int i=0; i<igraph_vector_size(&resultaux); i++){
		if((int)igraph_vector_e(&resultaux,i) <= 1){
			removed ++;
		}
	}

	/* inicializa el vector que contendra los nodos a remover */
	igraph_vector_init(&remaux, removed);

	int pos = 0; // permitira indicar en que posicion del array agregar un elemento
	/* se revisa los nodos de remover y se agregar al vector */
	for(int i=0; i<igraph_vector_size(&resultaux); i++){
		if((int)igraph_vector_e(&resultaux,i) <= 1){
			igraph_vector_set(&remaux, pos, i);
			pos++;
		}
	}	
	igraph_vector_set(&remaux, pos, -1); // convencion -1 que indique final de array

	/* proceso de eliminacion de los nodos seleccionados del grafo */
	igraph_vs_vector(&rem, &remaux); // crea el tipo de vector utilizado para llevar a cabo la remosion 

	igraph_delete_vertices(&graph, rem); // remueve los vertices y aristas asociadas

	igraph_vector_init(&result, 0);
	igraph_degree(&graph, &result, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS); // calcula los grados de los vertices que quedan en el grafo

	/*for(int i=0; i<igraph_vector_size(&result); i++){
		printf("Degree is  %10i, vertex %2i.\n", (int)igraph_vector_e(&result,i),i);
	}
	printf("Degree is  %10i, vertex %2i.\n", (int)igraph_vector_e(&result,393-16),393);*/

	/* calculo del nodo con mayor grado del 2-core */
	int max_node = igraph_vector_which_max(&result);
	printf("%i\n", max_node);

	/* remover nodo con mayor grado del 2-core */
	igraph_delete_vertices(&graph, igraph_vss_1(max_node));

	/* calculo de los nuevos grados de cada nodo */
	igraph_vector_init(&resultaux, 0);
	igraph_degree(&graph, &resultaux, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);

	printf("%li\n", igraph_vector_size(&resultaux));
	for(int i=0; i<igraph_vector_size(&resultaux); i++){
		printf("Degree is  %10i, vertex %2i.\n", (int)igraph_vector_e(&resultaux,i),i);
	}

}