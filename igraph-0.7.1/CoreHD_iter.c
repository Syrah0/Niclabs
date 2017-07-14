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

/* Se encarga de seleccionar los nodos que se deben eliminar para conformar el 2-core */
igraph_vector_t coreCal(igraph_vector_t *res){
	igraph_vector_t remaux;

	/* inicializa el vector que contendra los nodos a remover */
	igraph_vector_init(&remaux, 0);

	/* se revisa los nodos de remover y se agregar al vector */
	for(int i=0; i<igraph_vector_size(res); i++){
		if((int)igraph_vector_e(res,i) <= 1){
			igraph_vector_push_back(&remaux, i); // se agrega nodo que debe ser eliminado
		}
	}	
	return remaux;
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

/* Se encarga del proceso de eliminacion del nodo de mas grado y la actualizacion del 2-core y los grados del grafo */
igraph_t* deleteMaxNode(igraph_t *g){
	if(igraph_vcount(g) == 0){ // 2-core esta vacio
		return g;
	}
	else{ // eliminar nodo de mas grado y recalcular 2-core
		igraph_t* gaux;
		igraph_vector_t remaux, result;
		igraph_vs_t rem;
		int ret;

		igraph_vector_init(&result, 0);

		/* calcula los grados de los vertices que quedan en el grafo */
		igraph_degree(g, &result, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

		/* calculo del nodo con mayor grado del 2-core */
		int max_node = igraph_vector_which_max(&result);

		/* remover nodo con mayor grado del 2-core */
		ret = igraph_delete_vertices(g, igraph_vss_1(max_node));
		igraph_vector_destroy(&result);

		/* Calculo de la componente conexa mas grande */
		gaux = component(g);
		igraph_destroy(g);

		/* calculo de los nuevos grados de cada nodo */
		igraph_vector_init(&result, 0);
		igraph_degree(gaux, &result, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

		/* actualizar el 2-core */
		/* calcula que nodos se deben eliminar para obtener el 2-core */		
		remaux  = coreCal(&result);

		igraph_vector_destroy(&result);
		igraph_vs_vector(&rem, &remaux); // crea el tipo de vector utilizado para llevar a cabo la remosion 

		/* proceso de eliminacion de los nodos seleccionados del grafo */
		ret = igraph_delete_vertices(gaux, rem); // remueve los vertices y aristas asociadas -> se forma 2-core

		igraph_set_error_handler(igraph_error_handler_ignore);

		igraph_vector_destroy(&remaux);
		igraph_vs_destroy(&rem);

		/* Se verifica nuevamente el proceso anterior para el nuevo 2-core */
		deleteMaxNode(gaux);
	}
}

int main(){
	FILE *F;
	igraph_t graph;
	igraph_t* gaux;
	igraph_vector_t remaux, result;
	igraph_vs_t rem;
	int ret;

	F = fopen("red3.edges","r");
	igraph_read_graph_edgelist(&graph,F,0,0); // crea el grafo a partir del archivo con las conexiones
	fclose(F);

	igraph_simplify(&graph, 1, 1, 0); // se sacan posibles loops y aristas repetidas.

	/* Calculo de la componente conexa mas grande */
	gaux = component(&graph);
	igraph_destroy(&graph);

	/* calculo de los grados de cada nodo del grafo */
	igraph_vector_init(&result, 0);
	igraph_degree(gaux, &result, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS); 

	/* calcula que nodos se deben eliminar para obtener el 2-core */
	remaux  = coreCal(&result);
	igraph_vector_destroy(&result);

	igraph_vs_vector(&rem, &remaux); // crea el tipo de vector utilizado para llevar a cabo la remosion 

	/* proceso de eliminacion de los nodos seleccionados del grafo */
	ret = igraph_delete_vertices(gaux, rem); // remueve los vertices y aristas asociadas -> se forma 2-core

	igraph_set_error_handler(igraph_error_handler_ignore);

	igraph_vector_destroy(&remaux);
	igraph_vs_destroy(&rem);

	igraph_copy(&graph,deleteMaxNode(gaux)); // realizar paso 2 y 3 hasta que 2-core sea vacio

	fprintf(stderr, "%i\n", (int)igraph_vcount(&graph));
	igraph_destroy(gaux);

	printf("VACIO\n");
	return 0;
}