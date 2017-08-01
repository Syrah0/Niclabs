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

int main(){
	FILE *F, *G, *H;
	char filename[32];
	igraph_t graph, gaux, gaux2;
	igraph_vector_t degrees, nodes_aux;
	igraph_vs_t nodes_del;
	//double rest;
	double remove = 0.1; // multiplicador de porcentaje
	double rem_nodes = 0.0; // cantidad de nodos removidos
	int total_nodes; // cantidad de nodos del grafo original
	clock_t start, end, start_ini, end_ini;
	double time_used;
	double time_used_total = 0;

	G = fopen("Degree_times.csv","w"); // archivo que guardara los tiempo de ejecucion por iteracion
	H = fopen("Degree_iter.csv", "w"); // archivo que guardara los R-index (componente mas grande) por iteracion

	/* Se lee el archivo que contiene las conexiones de los nodos */
	F = fopen("red3.edges","r");
	igraph_read_graph_edgelist(&graph,F,0,0); // crea el grafo a partir del archivo con las conexiones
	fclose(F);
	total_nodes = igraph_vcount(&graph); // cantidad de nodos del grafo en analisis

//	start_ini = clock();
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
		//fprintf(stderr, "%s\n", output);
		fputs(output,G);
		putc(',',G);

		/* Proceso de escritura del nuevo grafo tras cierto porcentaje de eliminacion de nodos */
		rem_nodes += 1.0; // aumento cantidad de nodos removidos
		if(rem_nodes == ceil(total_nodes*remove)){
		    sprintf(filename,"grafo%d_Degree.edges",(total_nodes - (int)rem_nodes)); // nombre del archivo donde estaran los resultados
		    F = fopen(filename,"w");
		    igraph_write_graph_edgelist(&graph,F); // escritura
		    fclose(F);
		    remove += 0.1; // aumento porcentaje de eliminacion
		}
		fprintf(stderr, "%i\n", (int)rem_nodes);

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
//	end_ini = clock();

//	time_used = ((double) (end_ini - start_ini))/CLOCKS_PER_SEC;
	char output[50];		
		
	sprintf(output, "%f", time_used_total);
	fprintf(stderr, "%s\n", output);

	fclose(G);
	fclose(H);

	G = fopen("TiempoEjecución_Degree.txt","w");

	fputs(output,G);

	fclose(G);
	fprintf(stderr, "%i %i\n", igraph_vcount(&graph), igraph_ecount(&graph));

	//sprintf(filename,"grafoFinal_CI.edges"); // nombre del archivo donde estaran los resultados
	F = fopen("grafoFinal_Degree.edges","w");
        igraph_write_graph_edgelist(&graph,F); // escritura
	fclose(F);	

	printf("VACIO\n");

	return 0;
}
