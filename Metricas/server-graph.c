/* Funcion encargada de crear grafos aleatorios de a lo mas 1000 nodos */
/* Funcion paralelizable con threads por cada grafo que se cree */

#include <igraph.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
#include <pthread.h>
#include "kendall_sm.c"

typedef void *(*Thread_fun)(void *);

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
FILE *R;
float alpha;
int nodes, edges;
char *comp, *rad;

void createRandomGraph(char* filename, float alp, int nod, int edg){
	igraph_t graph;
	FILE *F;
	// generar power law graph
	igraph_static_power_law_game(&graph,nod,edg,alp,-1,1,0,0);
	F = fopen(filename,"w");
	igraph_write_graph_edgelist(&graph,F); // escritura

	fclose(F);

}

void* server(int pos){ 
	char filename[32];
	char **res;
	char *argv[4];
	igraph_t graph;
	igraph_integer_t diam;
	FILE *H;

	sprintf(filename, "archivos/Grafos/redTest%d.edges", pos);
	createRandomGraph(filename,alpha,nodes,edges);

	// Corrobora radio de CI
	H = fopen(filename,"r");
	igraph_read_graph_edgelist(&graph,H,0,0);
	igraph_diameter(&graph,&diam,0,0,0,IGRAPH_UNDIRECTED,1);

	if(atoi(rad)<0 || atoi(rad)>diam){
		fprintf(stderr, "Radio de CI fuera de los limites. 0 <= l <= %d\n",diam );
		exit(0);
	}

	fclose(H);
	igraph_destroy(&graph);


	argv[0] = filename; // -> filename
	argv[1] = rad; // debe ser cambiable segun cliente -> rad
	argv[2] = "order";
	argv[3] = comp; // debe ser cambiable segun cliente -> comp

	// llamado a kendall
	res = initKendall(argv,pos);
	sprintf(filename, "redTest%d.edges", pos);

	pthread_mutex_lock(&mutex);
	fputs(filename,R);
	putc(',',R);
	fputs(argv[3],R);
	putc(',',R);

	for(int i = 0; i < 12; i++){
		fputs(res[i],R);
		if( (i%2 == 1) && (i != 11)){
			putc(',',R);
		}
	}
	putc('\n',R);

	pthread_mutex_unlock(&mutex);
	return NULL;
}

int main(int argc, char **argv){

	if(argc == 7){
		alpha = atof(argv[2]);
		nodes = atoi(argv[3]);
		edges = atoi(argv[4]);
		comp = argv[5];
		rad = argv[6];
	}
	else{
		fprintf(stderr, "Use: threads alpha nodes edges elem-comp rad-CI\n");
		exit(1);
	}

	if(alpha < 2){
		fprintf(stderr, "Error. Alpha >= 2\n");
		exit(1);
	}

	long thr = 0;
	long pos = 0;
	pthread_t thread[atoi(argv[1])];

	R = fopen("archivos/Kendall_Correlation.csv","a");

	for(int i = 0; i < atoi(argv[1]); i++){ // cada vez que se conecta un cliente
		pthread_t pid;
		long myNum;
		while(1){
			myNum = pos;
			FILE *P;
			char filename[300];
			sprintf(filename, "archivos/Grafos/redTest%d.edges", (int)myNum);
			P = fopen(filename,"r");
			if(!P){
				break;
			}
			pos++;
		}

		if(pthread_create(&pid, NULL, (Thread_fun)server, (void *)myNum) != 0){
			fprintf(stderr, "No se puede crear nuevo thread de cliente\n");

		}

		thread[thr] = pid;

		fprintf(stderr, "thread lanzado\n");

		thr++;
		pos++;

	}

	for(int i = 0; i < atoi(argv[1]); i++){
		pthread_join(thread[i], NULL);
	}

	fclose(R);
	return 0;
}

