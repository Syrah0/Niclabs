/* 

1-core : todos los nodos
2-core : remover todos los nodos con un arco o menos. Remover arcos asociados a ese nodo.

Algoritmo CoreHD:
1.- Obtener el 2-core -> obtener grados de cada nodo del 2-core
2.- Obtener nodo con mayor grado
3.- Remover dicho nodo -> actualizar 2-core y grados de sus nodos -> verificar nvo 2-core: Vacío -> tree-breaking (dismantling)
							   															   No vacío -> Volver a paso 2	

*/

// POR REVISAR
	// VER BIEN SI SE DEBE HACER SOLO EN LA COMP MAS GRANDE //
	// O CONSIDERANDO EL GRAFO COMPLETO //
	// VER SI EL NODO MAX SE ELIMINA DEL 2-CORE Y DE AHI TRABAJAR 
	// O DEL GRAFO ORIGINAL
	// VER TEMA DE SI SACAR LOS LOOPS INICIALES O NO
	// VER SI LOS GRADOS DE LOS NODOS SE CALCULAN CON LOOPS O NO
	// Y VER LA COMP CONEXA MAS GRANDE 
	// VER IMPORTANCIA DEL 2-CORE CON NODOS DE GRADO 0 !! ES NECESARIO ELIMINARLOS

#include <igraph.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/* Se encarga de seleccionar los nodos que se deben eliminar para conformar el 2-core */
igraph_vector_t coreCal(igraph_vector_t *res, int val){
	igraph_vector_t remaux;

	/* inicializa el vector que contendra los nodos a remover */
	igraph_vector_init(&remaux, 0);

	/* se revisa los nodos de remover y se agregar al vector */
	for(int i=0; i<igraph_vector_size(res); i++){
		if((int)igraph_vector_e(res,i) != 0 && (int)igraph_vector_e(res,i) < val){
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

void print_vector(igraph_vector_t *v, FILE *f) {
  long int i;
  for (i=0; i<igraph_vector_size(v); i++) {
    fprintf(f, " %li", (long int) VECTOR(*v)[i]);
  }
  fprintf(f, "\n");
}

int main(){
	FILE *F;
	char filename[32];
	igraph_t graph, gaux;
	igraph_vector_t remaux, result, edges, alledges;
	igraph_es_t rem;
	igraph_vector_t res;
	int coreVal = 2;
	double remove = 0.1; // multiplicador de porcentaje
	double rem_nodes = 0.0; // cantidad de nodos removidos
	int total_nodes;

	F = fopen("red3.edges","r");
	igraph_read_graph_edgelist(&graph,F,0,0); // crea el grafo a partir del archivo con las conexiones
	fclose(F);

	total_nodes = igraph_vcount(&graph);
	fprintf(stderr, "%i\n", total_nodes);

	igraph_copy(&gaux, &graph); // gaux representara los 2-core formados	
	
	while(1){

		/* calculo de los grados de cada nodo del grafo */
		igraph_vector_init(&alledges,0);
		igraph_vector_init(&result, 0);
		igraph_degree(&gaux, &result, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS); 

		/* calcula que lados se deben eliminar para obtener el 2-core */
		remaux  = coreCal(&result, coreVal);
		//print_vector(&remaux, stdout);

		if(igraph_vector_size(&remaux) == 0){
			break;
		}
		

		for(int i=0; i<igraph_vector_size(&remaux); i++){
			igraph_vector_init(&edges,0);
			igraph_incident(&gaux,&edges,igraph_vector_e(&remaux,i),IGRAPH_ALL);
			for(int j=0; j<igraph_vector_size(&edges); j++){
				igraph_vector_push_back(&alledges, igraph_vector_e(&edges,j));
			}
		}

		igraph_es_vector(&rem, &alledges); // crea el tipo de vector utilizado para llevar a cabo la remosion 
		
		/* proceso de eliminacion de los lados seleccionados del grafo */
		igraph_delete_edges(&gaux, rem);

		igraph_vector_destroy(&result);
		igraph_es_destroy(&rem);
		igraph_vector_destroy(&remaux);
		igraph_vector_destroy(&alledges);
	}

	while(igraph_ecount(&gaux) != 0){ // realizar paso 2 y 3 hasta que 2-core sea vacio

		igraph_vector_init(&result, 0);
		igraph_degree(&gaux, &result, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS); // calcula los grados de los vertices que quedan en el grafo
		igraph_destroy(&gaux);

		/* calculo del nodo con mayor grado del 2-core */
		int max_node = igraph_vector_which_max(&result);

		//fprintf(stderr, "%i\n", max_node);
		/* remover nodo con mayor grado del grafo original */
		igraph_delete_vertices(&graph, igraph_vss_1(max_node));
		igraph_vector_destroy(&result);

		rem_nodes += 1.0;
		if(rem_nodes == ceil(total_nodes*remove)){
			sprintf(filename,"grafo%d_CoreHD.edges",(total_nodes - (int)rem_nodes));
		    F = fopen(filename,"w");
		    igraph_write_graph_edgelist(&graph,F);
		    fclose(F);
		    fprintf(stderr,"file: %s with %lf maxCC-nodes\n",filename,remove);
		    remove += 0.1;
		}
		fprintf(stderr, "%lf %lf\n", rem_nodes, ceil(total_nodes*(remove-0.1)));

		igraph_copy(&gaux,&graph); // para obtener el siguiente 2-core

		/* actualizar el 2-core */
		/* calcula que nodos se deben eliminar para obtener el 2-core */

		while(1){
			//fprintf(stderr, "ciclo2\n");
			/* calculo de los grados de cada nodo del grafo */
			igraph_vector_init(&result, 0);
			igraph_vector_init(&alledges,0);
			igraph_degree(&gaux, &result, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS); 
			
			/* calcula que lados se deben eliminar para obtener el 2-core */
			remaux  = coreCal(&result, coreVal);
			//print_vector(&remaux, stdout);

			if(igraph_vector_size(&remaux) == 0){
				break;
			}

			for(int i=0; i<igraph_vector_size(&remaux); i++){
				igraph_vector_init(&edges,0);
				igraph_incident(&gaux,&edges,igraph_vector_e(&remaux,i),IGRAPH_ALL);
				for(int j=0; j<igraph_vector_size(&edges); j++){
					igraph_vector_push_back(&alledges, igraph_vector_e(&edges,j));
				}
			}


			if(igraph_vector_size(&alledges) != 0){ // se eliminan los nodos
				igraph_es_vector(&rem, &alledges); // crea el tipo de vector utilizado para llevar a cabo la remosion 	
				/* proceso de eliminacion de los lado seleccionados del grafo */
				/* esto permite forma el 2-core sin eliminar los nodos y asi no borrar nodos erroneos del grafo original */
				igraph_delete_edges(&gaux, rem); // remueve los lados de los nodos que se deberian eliminar para el 2-core
				
			}
			igraph_es_destroy(&rem);
			igraph_vector_destroy(&result);
			igraph_vector_destroy(&remaux);
			igraph_vector_destroy(&alledges);
		}
		//print_vector(&alledges, stdout);
		//igraph_vector_destroy(&result);
	}

	fprintf(stderr, "%i %i\n", (int)igraph_vcount(&gaux), (int)igraph_ecount(&gaux));
	igraph_destroy(&gaux);
	fprintf(stderr, "%i %i\n", (int)igraph_vcount(&graph), (int)igraph_ecount(&graph));
	printf("VACIO\n");

	/* realizar tree-breaking o desmantelamiento */
	return 0;
}