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


/* Se encarga de calcular los nodos en el borde de la vecindad de radio l del nodo "node" */
igraph_vector_t nodesToDistance(igraph_t *g, int l, int node){
	igraph_vector_t nodes, nodesExc, neigh;

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
		for(int j = 0; j < igraph_vector_size(&nodes); j++){ // calcula los vecinos de cada nodo con distancia menor a l a "node"
			igraph_vector_init(&neigh, 0);
			igraph_neighbors(g, &neigh, igraph_vector_e(&nodes, j), IGRAPH_ALL); // calcula los vecinos del nodo
			for(int k = 0; k < igraph_vector_size(&nodesExc); k++){ // descarto nodos ya visitados anteriormente
				for(int m = 0; m < igraph_vector_size(&neigh); m++){
					if(igraph_vector_e(&neigh,m) == igraph_vector_e(&nodesExc, k)){ // si el nodo ya habia sido visitado previamente
						igraph_vector_remove(&neigh, m); // elimino dicho nodo
					}
				}
			}
			for(int k = 0; k < igraph_vector_size(&neigh); k++){ // agrego cada nodo en el borde de la vecindad calculada
				igraph_vector_push_back(&nodesExc, igraph_vector_e(&neigh, k)); // agrego a la lista de nodos ya visitados
				igraph_vector_push_back(&nodesaux, igraph_vector_e(&neigh, k)); // agrego a la lista de nodos pertenecientes al borde de la vecindad
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

int main(){
	FILE *F;
	char filename[32];
	igraph_t graph;
	igraph_vector_t degrees, nodes, CIvalues;
	int l = 2; // radio de la vecindad alegido
	double rest;
	double remove = 0.1; // multiplicador de porcentaje
	double rem_nodes = 0.0; // cantidad de nodos removidos
	int total_nodes; // cantidad de nodos del grafo original

	/* Se lee el archivo que contiene las conexiones de los nodos */
	F = fopen("red3.edges","r");
	igraph_read_graph_edgelist(&graph,F,0,0); // crea el grafo a partir del archivo con las conexiones
	fclose(F);

	total_nodes = igraph_vcount(&graph); 

	/* Calcular los grados de los nodos del grafo */
	igraph_vector_init(&degrees, 0);
	igraph_degree(&graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS); 

	double k = igraph_vector_sum(&degrees)/igraph_vector_size(&degrees); // grado promedio del grafo

	igraph_vector_destroy(&degrees);

	/* Calcular los grados de los nodos del grafo */
	while(1){
		/* Calcular los grados de los nodos del grafo */
		igraph_vector_init(&degrees, 0);
		igraph_degree(&graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS); 

		int N = igraph_vector_size(&degrees); // cantidad de nodos del grafo

		/* Calcular el CI de cada nodo */
		igraph_vector_init(&CIvalues, 0);
		for(int i = 0; i < N; i++){
			nodes = nodesToDistance(&graph, l, i); // calculo los nodos a distancia l del nodo "i"
			int CI = (igraph_vector_e(&degrees, i) - 1);
			int sum = 0;
			/* calculo de la sumatoria de la ecuacion CI a partir de los nodos en el borde de la vecindad de radio l */
			for(int j = 0; j < igraph_vector_size(&nodes); j++){ 
				sum += (igraph_vector_e(&degrees, igraph_vector_e(&nodes, j)) - 1);
			}
			CI *= sum; 
			igraph_vector_push_back(&CIvalues, CI); // agrega valor CI calculado al vector que contendra todos los CI de cada nodo
		}

		double y = (1.0/(l+1)); // calcula el exponente de la condicion de termino del algoritmo
		double x = (igraph_vector_sum(&CIvalues)/(N * k)); // calcula la base de la condicion de termino del algoritmo
		rest = pow(x,y); // calcula la condicion de termino

		if(rest <= 1){ // verifica condicion de termino del algoritmo
			break;
		}

		int max_node = igraph_vector_which_max(&CIvalues); // obtiene el nodo con mayor CI
		igraph_delete_vertices(&graph, igraph_vss_1(max_node)); // eliminacion del nodo con mayor CI

		/* Proceso de escritura del nuevo grafo tras cierto porcentaje de eliminacion de nodos */
		rem_nodes += 1.0; // aumento cantidad de nodos removidos
		if(rem_nodes == ceil(total_nodes*remove)){
			sprintf(filename,"grafo%d_CI2.edges",(total_nodes - (int)rem_nodes)); // nombre del archivo donde estaran los resultados
		    F = fopen(filename,"w");
		    igraph_write_graph_edgelist(&graph,F); // escritura
		    fclose(F);
		    remove += 0.1; // aumento porcentaje de eliminacion
		}
		fprintf(stderr, "%i\n", (int)rem_nodes);
	}
	fprintf(stderr, "%i %i\n", igraph_vcount(&graph), igraph_ecount(&graph));

	printf("VACIO\n");

	return 0;
}