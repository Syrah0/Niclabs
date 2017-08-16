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

int init_Betweenness(const char * name){
	FILE *F, *G, *H;
	char filename[32];
	igraph_t graph;
	igraph_vector_t betweenness, degrees, nodes;
	double remove = 0.1; // multiplicador de porcentaje
	double rem_nodes = 0.0; // cantidad de nodos removidos
	int total_nodes; // cantidad de nodos del grafo original
	clock_t start, end, start_ini, end_ini;
	double time_used;
	double time_used_total = 0;

	G = fopen("Betweenness_times.csv","w"); // archivo que guardara los tiempo de ejecucion por iteracion
	H = fopen("Betweenness_iter.csv", "w"); // archivo que guardara los R-index (componente mas grande) por iteracion

	/* Se lee el archivo que contiene las conexiones de los nodos */
	F = fopen(name,"r");
	if(F == NULL){
		fprintf(stderr, "Error al abrir el archivo %s\n", name);
		exit(1);
	}
	igraph_read_graph_edgelist(&graph,F,0,0); // crea el grafo a partir del archivo con las conexiones
	fclose(F);
	
	igraph_simplify(&graph,1,0,/*edge_comb=*/ 0);

	igraph_vector_init(&nodes,0);
	total_nodes = igraph_vcount(&graph); // cantidad de nodos del grafo en analisis
	int del_nodes[total_nodes];
	for(int i = 0; i < total_nodes; i++){
		del_nodes[i] = i;
	}

	int node, giant_comp;
	int iter = 1;
	while(1){
		start = clock();
		igraph_vector_init(&betweenness, 0);
		igraph_betweenness(&graph, &betweenness, igraph_vss_all(), 0,/*weights=*/ 0, 0);  // calculo del betweenness de cada nodo del grafo
		node = igraph_vector_which_max(&betweenness);

		//print_vector(&betweenness,stdout);

		igraph_delete_vertices(&graph,igraph_vss_1(node));
		igraph_vector_destroy(&betweenness);

		/* agregar nodo removido a la lista */
		int rest = 0;
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

		end = clock();
		giant_comp = max_component(&graph);

		time_used = ((double) (end - start))/CLOCKS_PER_SEC;
		time_used_total += time_used;
		char output[50];		

		sprintf(output, "%f", time_used);
		fprintf(stderr, "%s\n", output);
		fputs(output,G);
		putc(',',G);

		/* Proceso de escritura del nuevo grafo tras cierto porcentaje de eliminacion de nodos */
		rem_nodes += 1.0; // aumento cantidad de nodos removidos
		if(rem_nodes == ceil(total_nodes*remove)){
		    sprintf(filename,"grafo%d_Betweenness.edges",(total_nodes - (int)rem_nodes)); // nombre del archivo donde estaran los resultados
		    F = fopen(filename,"w");
		    igraph_write_graph_edgelist(&graph,F); // escritura
		    fclose(F);
		    remove += 0.1; // aumento porcentaje de eliminacion
		}

		sprintf(output,"%d", (int)iter);
		fputs(output,H);

		putc(',',H);

		sprintf(output, "%d", giant_comp);

		fputs(output,G);
		putc('\n',G);

		fputs(output,H);
		putc('\n',H);

		if(giant_comp == 1){
			break;
		}
		iter++;		
	}
	//print_vector(&nodes,stdout);
	char output[50];	

	F = fopen("removedNodes_Betweenness.txt","w");
	for(int i = 0; i < igraph_vector_size(&nodes)-1; i++){
		sprintf(output, "%d\n", (int)igraph_vector_e(&nodes,i));
		fputs(output,F);
	}
	sprintf(output, "%d", (int)igraph_vector_e(&nodes,igraph_vector_size(&nodes)-1));
	fputs(output,F);
	fclose(F);	
		
	sprintf(output, "%f", time_used_total);
	fprintf(stderr, "%s\n", output);

	fclose(G);
	fclose(H);

	G = fopen("TiempoEjecución_Betweenness.txt","w");

	fputs(output,G);

	fclose(G);
	fprintf(stderr, "%i %i\n", igraph_vcount(&graph), igraph_ecount(&graph));

	F = fopen("grafoFinal_Betweenness.edges","w");
    igraph_write_graph_edgelist(&graph,F); // escritura
	fclose(F);	

	printf("VACIO\n");

	return 0;
}

int main(int argc, const char * argv[]){

	if(argc != 2){
		fprintf(stderr, "Por favor ingresar el nombre del archivo que contiene el grafo con su extensión\n");
		exit(1);
	}

	init_Betweenness(argv[1]);

	return 0;
}
