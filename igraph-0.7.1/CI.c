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

/* se encarga de calcular los nodos dentro de la vecindad de radio l del nodo "node" */
igraph_vector_t neighborhood(igraph_t *g, int l, int node){
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
		/* por cada nodo perteneciente al borde la de vecindad anterior se calculan sus vecinos
		esto permite ampliar el tamaño de la vecindad */
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
				igraph_vector_push_back(&nodesExc, igraph_vector_e(&neigh, k)); // agrego a la lista de nodos ya visitados -> pertenecen a la vecindad de radio l
				igraph_vector_push_back(&nodesaux, igraph_vector_e(&neigh, k)); // agrego a la lista de nodos pertenecientes al borde de la vecindad
			}
			igraph_vector_destroy(&neigh);
		}

		igraph_vector_destroy(&nodes); // destruye vector para realizar copia
		igraph_vector_copy(&nodes, &nodesaux);

		igraph_vector_destroy(&nodesaux);
		igraph_vector_destroy(&neigh);
	}

	igraph_vector_destroy(&nodes);
	igraph_vector_remove(&nodesExc,0); // remuevo el nodo central "node"
	return nodesExc;
}

/* Se encarga de calcular el tamaño de la componente conexa mas grande */
int max_component(igraph_t *g){
	igraph_vector_ptr_t complist;

	/* Se inicializa vector que contendra todas las componentes conexas del grafo */
	igraph_vector_ptr_init(&complist, 0);
	igraph_decompose(g,&complist,IGRAPH_WEAK,-1,2); // se realiza la descomposicion del grafo en sus componentes

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
	igraph_t graph, gaux;
	igraph_vector_t degrees, nodes, CIvalues, nodes_neigh, edges;
	igraph_es_t rem;
	int l = 2; // radio de la vecindad elegido
	double rest;
	double remove = 0.1; // multiplicador de porcentaje
	double rem_nodes = 0.0; // cantidad de nodos removidos
	int total_nodes; // cantidad de nodos del grafo original
	clock_t start, end;
	double time_used;

	G = fopen("CI_times.csv","w"); // archivo que guardara los tiempo de ejecucion por iteracion
	H = fopen("CI_iter.csv", "w"); // archivo que guardara los R-index (componente mas grande) por iteracion

	/* Se lee el archivo que contiene las conexiones de los nodos */
	F = fopen("red3.edges","r");
	igraph_read_graph_edgelist(&graph,F,0,0); // crea el grafo a partir del archivo con las conexiones
	fclose(F);

	start = clock();
	/* Calcular los grados de los nodos del grafo */
	igraph_vector_init(&degrees, 0);
	igraph_degree(&graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);  // calculo de los grados de cada nodo del grafo

	double k = igraph_vector_sum(&degrees)/igraph_vector_size(&degrees); // grado promedio del grafo

	int N = igraph_vcount(&graph); // cantidad de nodos del grafo en analisis
	total_nodes = N;

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

	while(rest > 1){
		int max_node = igraph_vector_which_max(&CIvalues); // obtiene el nodo con mayor CI

		/* Calcular el CI de cada nodo */
		nodes_neigh = neighborhood(&graph, l+1, max_node); // veo nodos dentro de vecindad l+1 del nodo a eliminar
		
		igraph_vector_init(&edges,0);
	    igraph_incident(&graph,&edges,max_node,IGRAPH_ALL); // calcula los lados incidentes al nodo a eliminar para aislarlo
			
		igraph_es_vector(&rem, &edges); 
		igraph_delete_edges(&graph, rem); // eliminacion de los vertices incidentes al nodo a eliminar

		igraph_vector_destroy(&edges);
		igraph_es_destroy(&rem);
		igraph_vector_destroy(&degrees);

		// HACER COPIA A GAUX PARA BORRAR NODO AHI, LUEGO COPIAR GAUX A GRAPH 
		igraph_copy(&gaux, &graph);
		igraph_delete_vertices(&gaux, igraph_vss_1(max_node)); // eliminacion del nodo con mayor CI

		end = clock();
		time_used = ((double) (end - start))/CLOCKS_PER_SEC;
		char output[50];		
		
		sprintf(output, "%f", time_used);
		//fprintf(stderr, "%s\n", output);

		fputs(output,G);
		putc(',',G);

		/* Proceso de escritura del nuevo grafo tras cierto porcentaje de eliminacion de nodos */
		rem_nodes += 1.0; // aumento cantidad de nodos removidos
		if(rem_nodes == ceil(total_nodes*remove)){
			sprintf(filename,"grafo%d_CI.edges",(total_nodes - (int)rem_nodes)); // nombre del archivo donde estaran los resultados
		    F = fopen(filename,"w");
		    igraph_write_graph_edgelist(&gaux,F); // escritura
		    fclose(F);
		    remove += 0.1; // aumento porcentaje de eliminacion
		}
		fprintf(stderr, "%i\n", (int)rem_nodes);

		/* fin iteracion */
		sprintf(output,"%d", (int)rem_nodes);
		fputs(output,H);

		putc(',',H);

		int giant_comp = max_component(&gaux);
		sprintf(output, "%d", giant_comp);

		fputs(output,G);
		putc('\n',G);

		fputs(output,H);
		putc('\n',H);

		/* inicio siguiente iteracion */

		start = clock();

		/* Calcular los grados de los nodos del grafo */
		igraph_vector_init(&degrees, 0);
		igraph_degree(&graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS); 

		N = igraph_vcount(&graph) - 1; // cantidad de nodos del grafo (-1 porque no se debe considerar el nodo que se desea eliminar)

		/* actualizacion de los valores CI solo de los nodos dentro de la vecindad de radio l+1 del nodo a eliminar */
		for(int i = 0; i < igraph_vector_size(&nodes_neigh); i++){
			nodes = nodesToDistance(&graph, l, igraph_vector_e(&nodes_neigh,i)); // calcula los nodos en el borde de la vecindad de radio l
			int CI = (igraph_vector_e(&degrees, igraph_vector_e(&nodes_neigh,i)) - 1);
			int sum = 0;
			/* calculo de la sumatoria de la ecuacion CI a partir de los nodos en el borde de la vecindad de radio l */
			for(int j = 0; j < igraph_vector_size(&nodes); j++){
				sum += (igraph_vector_e(&degrees, igraph_vector_e(&nodes, j)) - 1);
			}
			CI *= sum;
			igraph_vector_set(&CIvalues,igraph_vector_e(&nodes_neigh,i),CI); // actualiza el valor de CI del nodo i dentro del vector
		}
		
		igraph_vector_remove(&CIvalues, max_node); // elimina el valor de CI del nodo a eliminar
		igraph_copy(&graph,&gaux);
		igraph_destroy(&gaux);

		//igraph_delete_vertices(&graph, igraph_vss_1(max_node)); // eliminacion del nodo con mayor CI

		/* Proceso de escritura del nuevo grafo tras cierto porcentaje de eliminacion de nodos */
/*		rem_nodes += 1.0; // aumento cantidad de nodos removidos
		if(rem_nodes == ceil(total_nodes*remove)){
			sprintf(filename,"grafo%d_CI.edges",(total_nodes - (int)rem_nodes)); // nombre del archivo donde estaran los resultados
		    F = fopen(filename,"w");
		    igraph_write_graph_edgelist(&graph,F); // escritura
		    fclose(F);
		    remove += 0.1; // aumento porcentaje de eliminacion
		}
		fprintf(stderr, "%i\n", (int)rem_nodes);

		/* calculo de la nueva condicion de termino */
		x = (igraph_vector_sum(&CIvalues)/(N * k)); // calculo de la nueva base
		rest = pow(x,y);
	}
	fclose(H);
	fclose(G);
	fprintf(stderr, "%i %i\n", igraph_vcount(&graph), igraph_ecount(&graph));
	printf("VACIO\n");

	return 0;
}