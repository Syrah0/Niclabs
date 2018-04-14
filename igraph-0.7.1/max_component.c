#include <igraph.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Se encarga de calcular el tamaño de la componente conexa mas grande */
int max_component(igraph_t *g){
	igraph_vector_ptr_t complist;

	/* Se inicializa vector que contendra todas las componentes conexas del grafo */
	igraph_vector_ptr_init(&complist, 0);
	igraph_decompose(g,&complist,IGRAPH_WEAK,-1,0); // se realiza la descomposicion del grafo en sus componentes

	/* Se verifica cual es la componente mas grande */
	int max = 0;
	for(int i=1; i<igraph_vector_ptr_size(&complist); i++){
		if(igraph_vcount(VECTOR(complist)[max]) < igraph_vcount(VECTOR(complist)[i])){ // Hay una mas grande que la componente actual
			max = i;
		}
	}
	return (int)(igraph_vcount(VECTOR(complist)[max])); // componente conexa mas grande
}
