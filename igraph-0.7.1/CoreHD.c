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
	igraph_t graph;
	igraph_vector_t result;
	F = fopen("red3.edges","r");
	igraph_read_graph_edgelist(&graph,F,0,0);
	fclose(F);

	igraph_vector_init(&result, 0);
	
	igraph_degree(&graph, &result, igraph_vss_all(), IGRAPH_ALL,
	               IGRAPH_LOOPS);

	for(int i=0; i<igraph_vector_size(&result); i++){
		printf("Degree is  %10i, vertex %2i.\n", (int)igraph_vector_e(&result,i),i);
	}
}