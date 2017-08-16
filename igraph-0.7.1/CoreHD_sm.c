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
#include <math.h>
#include <time.h>
#include "max_component.c"


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

int init_CoreHD(const char *name){
	FILE *F, *G, *H;
	char filename[32];
	igraph_t graph, gaux;
	igraph_vector_t remaux, result, edges, alledges, res, nodes;
	igraph_es_t rem;
	int coreVal = 2; // k-core deseado
	double remove = 0.1; // multiplicador de porcentaje
	double rem_nodes = 0.0; // cantidad de nodos removidos
	int total_nodes; // cantidad de nodos del grafo original
	clock_t start, end, start_ini, end_ini;
	double time_used;
	double time_used_total = 0;

	G = fopen("CoreHD_times.csv","w"); // archivo que guardara los tiempos de ejecucion por iteracion
	H = fopen("CoreHD_iter.csv", "w"); // archivo que guardara los R-index (componente mas grande) por iteracion

	start = clock();
	/* se lee el archivo con los datos del grafo */
	F = fopen(name,"r");
	if(F == NULL){
		fprintf(stderr, "Error al abrir el archivo %s\n", name);
		exit(1);
	}
	igraph_read_graph_edgelist(&graph,F,0,0); // crea el grafo a partir del archivo con las conexiones
	fclose(F);

	igraph_simplify(&graph,1,0,/*edge_comb=*/ 0);

	total_nodes = igraph_vcount(&graph);
	
	igraph_vector_init(&nodes,0);
	int del_nodes[total_nodes];
	for(int i = 0; i < total_nodes; i++){
		del_nodes[i] = i;
	}

	//fprintf(stderr, "%i\n", total_nodes);

	igraph_copy(&gaux, &graph); // gaux representara los 2-core formados	
	
	while(1){

		/* calculo de los grados de cada nodo del grafo */
		igraph_vector_init(&alledges,0);
		igraph_vector_init(&result, 0);
		igraph_degree(&gaux, &result, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS); 

		/* calcula que lados se deben eliminar para obtener el 2-core */
		remaux  = coreCal(&result, coreVal);

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
		//fprintf(stderr, "max: %d\n", max_node);

		/* remover nodo con mayor grado del grafo original */
		igraph_delete_vertices(&graph, igraph_vss_1(max_node));

		/* agregar nodo removido a la lista */
		int rest = 0;
		for(int i = 0; i < total_nodes; i++){
			if(del_nodes[i] == max_node){
				// agrego a lista
				igraph_vector_push_back(&nodes,i);
				del_nodes[i] = -1;
				rest = 1;
			}
			else{
				del_nodes[i] -= rest;
			}
		}

		igraph_vector_destroy(&result);

		end = clock();
		time_used = ((double) (end - start))/CLOCKS_PER_SEC;
		time_used_total += time_used;
		char output[50];		
		
		sprintf(output, "%f", time_used);
		//fprintf(stderr, "%s\n", output);

		fputs(output,G);
		putc(',',G);

		/* Proceso de escritura del nuevo grafo tras cierto porcentaje de eliminacion de nodos */
		rem_nodes += 1.0; // aumento cantidad de nodos removidos
		if(rem_nodes == ceil(total_nodes*remove)){
			sprintf(filename,"grafo%d_CoreHD.edges",(total_nodes - (int)rem_nodes)); // nombre del archivo donde estaran los resultados
		    F = fopen(filename,"w");
		    igraph_write_graph_edgelist(&graph,F); // escritura
		    fclose(F);
		    remove += 0.1; // aumento porcentaje de eliminacion
		}

		/* fin iteracion */
		sprintf(output,"%d", (int)rem_nodes);
		fputs(output,H);

		putc(',',H);

		int giant_comp = max_component(&graph);
		sprintf(output, "%d", giant_comp);

		fputs(output,G);
		putc('\n',G);

		fputs(output,H);
		putc('\n',H);

		/* inicio siguiente iteracion */
		start = clock();
		igraph_copy(&gaux,&graph); // para obtener el siguiente 2-core

		/* actualizar el 2-core */
		/* calcula que nodos se deben eliminar para obtener el 2-core */

		while(1){
			/* calculo de los grados de cada nodo del grafo */
			igraph_vector_init(&result, 0);
			igraph_vector_init(&alledges,0);
			igraph_degree(&gaux, &result, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS); 
			
			/* calcula que lados se deben eliminar para obtener el 2-core */
			remaux  = coreCal(&result, coreVal);

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
	}
	char output[50];		
		
	sprintf(output, "%f", time_used_total);
	//fprintf(stderr, "%s\n", output);

	fclose(G);
	fclose(H);

	G = fopen("TiempoEjecución_CoreHD.txt","w");

	fputs(output,G);

	fclose(G);	

	igraph_destroy(&gaux);
	//fprintf(stderr, "%i %i\n", (int)igraph_vcount(&graph), (int)igraph_ecount(&graph));
	
	F = fopen("grafoFinal_CoreHD.edges","w");
    igraph_write_graph_edgelist(&graph,F); // escritura
	fclose(F);

	/* realizar tree-breaking o desmantelamiento */
	fprintf(stderr, "Desmantelamiento\n");

	while(1){
		int giant_comp = max_component(&graph);
		if(giant_comp == 1){
			break;
		}

		igraph_vector_init(&result, 0);
		igraph_degree(&graph, &result, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS); 

		int max_node = igraph_vector_which_max(&result);
		igraph_delete_vertices(&graph, igraph_vss_1(max_node));
		int rest = 0;
		for(int i = 0; i < total_nodes; i++){
			if(del_nodes[i] == max_node){
				// agrego a lista
				igraph_vector_push_back(&nodes,i);
				del_nodes[i] = -1;
				rest = 1;
			}
			else{
				del_nodes[i] -= rest;
			}
		}
	}	

	F = fopen("removedNodes_CoreHD.txt","w");
	for(int i = 0; i < igraph_vector_size(&nodes)-1; i++){
		sprintf(output, "%d\n", (int)igraph_vector_e(&nodes,i));
		fputs(output,F);
	}
	sprintf(output, "%d", (int)igraph_vector_e(&nodes,igraph_vector_size(&nodes)-1));
	fputs(output,F);
	fclose(F);

	fprintf(stderr, "%i %i\n", (int)igraph_vcount(&graph), (int)igraph_ecount(&graph));
	printf("VACIO\n");

	return 0;
}
