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
#include <math.h>
#include <time.h>
//#include "max_component.c"

int init_Degree(char * name, int nodeComp, int pos){
	FILE *F, *G, *H;
	char filename[300];
	igraph_t graph, gaux, gaux2;
	igraph_vector_t degrees, nodes_aux, nodes;
	igraph_vs_t nodes_del;
	double remove = 0.1; // multiplicador de porcentaje
	double rem_nodes = 0.0; // cantidad de nodos removidos
	int total_nodes; // cantidad de nodos del grafo original
	clock_t start, end, start_ini, end_ini;
	double time_used;
	double time_used_total = 0;

//	fprintf(stderr, "Corriendo Degree\n");

/* GENERACION DE ARCHIVOS
	sprintf(filename, "Degree_times_graph%d.csv",pos);
	G = fopen(filename,"w"); // archivo que guardara los tiempo de ejecucion por iteracion

	sprintf(filename, "Degree_iter_graph%d.csv",pos);
	H = fopen(filename, "w"); // archivo que guardara los R-index (componente mas grande) por iteracion
*/
	/* Se lee el archivo que contiene las conexiones de los nodos */
	F = fopen(name,"r");
	if(F == NULL){
		fprintf(stderr, "Error al abrir el archivo %s\n", name);
		exit(1);
	}
	igraph_read_graph_edgelist(&graph,F,0,0); // crea el grafo a partir del archivo con las conexiones
	fclose(F);
	igraph_simplify(&graph,1,0,/*edge_comb=*/ 0);
	total_nodes = igraph_vcount(&graph); // cantidad de nodos del grafo en analisis

	igraph_vector_init(&nodes,0);
	int del_nodes[total_nodes];
	for(int i = 0; i < total_nodes; i++){
		del_nodes[i] = i;
	}

	start = clock();
	int node, giant_comp, degree;
	int iter = 1;
	int values[total_nodes];

	igraph_vector_init(&nodes_aux, 0);
	igraph_vector_init(&degrees, 0);
	igraph_degree(&graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);  // calculo de los grados de cada nodo del grafo
	
	for(int j=0; j < total_nodes; j++){
		degree = (int)igraph_vector_max(&degrees);
		node = (int)igraph_vector_which_max(&degrees);
		values[j] = node;
		igraph_vector_set(&degrees,node,-1);
	}
	igraph_vector_destroy(&degrees);
	igraph_copy(&gaux,&graph);
	igraph_copy(&gaux2,&graph);

	while(1){

		if(iter != 1) start = clock();

		node = values[iter-1]; // sera el mayor de cada iteracion
		igraph_vector_push_back(&nodes_aux,node);
		igraph_vs_vector(&nodes_del,&nodes_aux);
		igraph_delete_vertices(&gaux,nodes_del);

		igraph_copy(&graph,&gaux);
		igraph_copy(&gaux,&gaux2);

		end = clock();
		giant_comp = max_component(&graph);

		time_used = ((double) (end - start))/CLOCKS_PER_SEC;
		time_used_total += time_used;
		char output[50];		

		sprintf(output, "%f", time_used);
		
/* DESCOMENTAR SI SE USA G
		fputs(output,G);
		putc(',',G);
*/
		/* Proceso de escritura del nuevo grafo tras cierto porcentaje de eliminacion de nodos */
		rem_nodes += 1.0; // aumento cantidad de nodos removidos
		
/* GENERACION ARCHIVOS POR PORCENTAJE REMOVIDO
		if(rem_nodes == ceil(total_nodes*remove)){
		    sprintf(filename,"grafo%d_Degree_graph%d.edges",(total_nodes - (int)rem_nodes),pos); // nombre del archivo donde estaran los resultados
		    F = fopen(filename,"w");
		    igraph_write_graph_edgelist(&graph,F); // escritura
		    fclose(F);
		    remove += 0.1; // aumento porcentaje de eliminacion
		}
*/
		//fprintf(stderr, "%i\n", (int)rem_nodes);
		
/* DESCOMENTAR SI SE USA G Y/O H
		sprintf(output,"%d", (int)iter);
		fputs(output,H);

		putc(',',H);

		sprintf(output, "%d", giant_comp);

		fputs(output,G);
		putc('\n',G);

		fputs(output,H);
		putc('\n',H);
*/
		if(giant_comp == 1 || (int)rem_nodes == nodeComp){
			break;
		}
		iter++;		
	}
	char output[50];	

	for(int i = 0; i < igraph_vector_size(&nodes_aux); i++){
		int rest = 0;
		for(int j = 0; j < total_nodes; j++){
			if(j == (int)igraph_vector_e(&nodes_aux,i)){
				// agrego a lista
				igraph_vector_push_back(&nodes,j);
				del_nodes[j] = -1;
				rest = 1;
			}
			else{
				if(del_nodes[j] != -1){
					del_nodes[j] -= rest;
				}
			}
		}
	}

//	fprintf(stderr, "%i %i\n", (int)igraph_vcount(&graph), (int)igraph_ecount(&graph));
	
//	fprintf(stderr, "Eliminación nodos grado 0\n");	

	/* Eliminacion nodos de grado cero */
	while(igraph_vcount(&graph) > 0 && (int)rem_nodes < nodeComp){
		igraph_vector_init(&degrees,0);
		igraph_degree(&graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS); 
		node = igraph_vector_which_max(&degrees);

		igraph_delete_vertices(&graph,igraph_vss_1(node));
		igraph_vector_destroy(&degrees);

		/* agregar nodo removido a la lista */
		int rest = 0;
		rem_nodes += 1.0;
		for(int i = 0; i < total_nodes; i++){
			if(del_nodes[i] == node){
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

	if(pos == -1){
		F = fopen("removedNodes_Degree.txt","w");
	}
	else{
		sprintf(filename, "archivos/Degree/removedNodes_Degree_graph%d.txt",pos);
		F = fopen(filename,"w");
	}
	
	for(int i = 0; i < igraph_vector_size(&nodes)-1; i++){
		sprintf(output, "%d\n", (int)igraph_vector_e(&nodes,i));
		fputs(output,F);
	}
	sprintf(output, "%d", (int)igraph_vector_e(&nodes,igraph_vector_size(&nodes)-1));
	fputs(output,F);
	fclose(F);
		
	//sprintf(output, "%f", time_used_total);
	//fprintf(stderr, "%s\n", output);
	
/* DESCOMENTAR SI SE USA G Y/O H
	fclose(G);
	fclose(H);

/* GENERACION ARCHIVO DE TIEMPO DE EJECUCION
	sprintf(filename, "TiempoEjecución_Degree_graph%d.txt",pos);

	G = fopen(filename,"w");

	fputs(output,G);

	fclose(G);
*/

	//fprintf(stderr, "%i %i\n", igraph_vcount(&graph), igraph_ecount(&graph));

/*	GENERACION ARCHIVO GRAFO FINAL

	sprintf(filename, "grafoFinal_Degree_graph%d.edges",pos);
	F = fopen(filename,"w");
    igraph_write_graph_edgelist(&graph,F); // escritura
	fclose(F);	
*/
	//printf("VACIO\n");

	return 0;
}
