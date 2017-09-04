// KENDALL
// Xi are the node from one metric and Yi are the node from the other metrica
// Xi and Yi are differents between them

// Calculate the number of concordant pairs and discordant pairs
/* Concordant:
	if Xi < Xj and Yi < Yj => (Xi, Yi) and (Xj, Yj) concordants pairs
	if Xi > Xj and Yi > Yj => (Xi, Yi) and (Xj, Yj) concordants pairs
	if Xi == Xj or Yi == Yj => any classification
	in other case, they are discordant

	=> r = (concordant - discordant)/ [n(n-1)/2] 
	=> -1 <= r <= 1 ... if X and Y are independent => r ~ 0

COMBINATION BETWEEN ALL PAIRS -> X1,Y1 ~ X2,Y2 - X3,Y3 - ... - XN,YN

=> N = TOTAL OF PAIRS
*/

// x1,y1 - x2,y2 - x3,y3 - x4,y4 - x5,y5 - x6,y6
//   0       1       2       3       4       5
// len = 6  => len-1 = 5
// i = 0 => j = 1 - 2 - 3 - 4 - 5
// i = 1 => j = 2 - 3 - 4 - 5
// i = 2 => j = 3 - 4 - 5 
// i = 3 => j = 4 - 5
// i = 4 => j = 5

#include <igraph.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "CoreHD_sm.c"
#include "degree_sm.c"
#include "Betweenness_sm.c"
#include "CI_sm.c"


double kendallCoef(int* X, int* Y, int len){

	int concordant = 0, discordant = 0;
	double res,div;
	for(int i = 0; i<(len-1); i++){ // i = 0 ... 4
		for(int j = (i+1); j<len; j++){ // j = i+1 ... 5
			if((X[i] > X[j] && Y[i] > Y[j]) || (X[i] < X[j] && Y[i] < Y[j])){
				concordant++;
			}
			else if(X[i] == X[j] || Y[i] == Y[j]){
				continue;
			}
			else{
				discordant++;
			}
		}
	}
	res = (concordant - discordant);
	div = (len*(len-1))/2;
	fprintf(stderr, "C=%d D=%d div=%lf\n", concordant,discordant,div);
	return res/div;
}

int CalculateKendallMetric(FILE *F, FILE *G, FILE *R, char *compare, int comp){
	igraph_vector_t nodes_met1, nodes_met2;
	int node, lenMet1, lenMet2, len, mode;
	char output[20];

	if(F == NULL || G == NULL){
		fprintf(stderr, "Error al abrir uno o ambos archivos\n");
		exit(1);
	}

	mode = fseek(F,0,SEEK_SET);
	if(mode){
		fprintf(stderr, "error fseek primer archivo\n");
	}

	mode = fseek(G,0,SEEK_SET);
	if(mode){
		fprintf(stderr, "error fseek segundo archivo\n");
	}

	igraph_vector_init(&nodes_met1,0);
	igraph_vector_init(&nodes_met2,0);

	// Copiar nodos de la metrica 1 en array
	while(fscanf(F, "%d", &node) > 0){
		igraph_vector_push_back(&nodes_met1,node);
	}

	// Copiar nodos de la metrica 2 en array
	while(fscanf(G, "%d", &node) > 0){
		igraph_vector_push_back(&nodes_met2,node);
	}

	lenMet1 = igraph_vector_size(&nodes_met1);
	lenMet2 = igraph_vector_size(&nodes_met2);

	//fprintf(stderr, "%d, %d\n", lenMet1,lenMet2);

	// revisar que largo de metrica es menor...
	// para hacer los pares segun dicha metrica
	if(comp < lenMet1 && comp < lenMet2){
		len = comp;
	}
	else if(lenMet1 < lenMet2){
		len = lenMet1;
	}
	else{
		len = lenMet2;
	}

	int met1[len];
	int met2[len];

	for(int i = 0; i < len; i++){
		met1[i] = (int)igraph_vector_e(&nodes_met1,i);
	}

	for(int i = 0; i < len; i++){
		met2[i] = (int)igraph_vector_e(&nodes_met2,i);
	}

	double kendall = kendallCoef(met1,met2,len);
	sprintf(output,"%f",kendall);

	fputs(compare,R);
	fputs(output,R);

	fprintf(stderr, "%lf\n", kendall);
}

int main(int argc, const char * argv[]){

	// VER INGRESO DE ARCHIVO CON EL GRAFO -> CALCULAR LOS NODOS REMOVIDOS 
	// CALCULAR LAS 4 METRICAS CON EL GRAFO (SIN MS CON CI)->  CREO QUE SON 6 COMPARACIONES -> 6 RESULTADOS (pero ver si se hace con CI o no!!!!)
	// VER COMO LLAMAR A LAS FUNCIONES YA HECHAS CON LAS METRICAS DESDE ACA -> VER SI SE PUEDE CON MAIN SINO CREAR FUN RUN LA CUAL ES LLAMADA DENTRO DE MAIN
	// HACER QUE EL MAIN DE LOS CODIGOS RECIBA EL NOMBRE DEL ARCHIVO CON EL GRAFO --> HACER ESTO
	// LISTO // CALCULAR LOS NODOS REMOVIDOS DE CADA ALGORITMO MENOS MINSUM 
	// LISTO // ARREGLAR COREHD CON PROCESO DE DESMANTELAMIENTO

	FILE *F, *G, *H, *I, *R;
	int rad = 0, comp = 0;
	const char* num = argv[2];
	const char* elem = argv[3];

	if(argc != 4){
		fprintf(stderr, "Ingresar: nombre del archivo que contiene el grafo, la vecindad de CI y cantidad de elementos a comparar\n");
		exit(1);
	}

	for(int i = 0; i < strlen(num); i++){
		rad = rad*10 + (num[i] - '0');
	}

	for(int i = 0; i < strlen(elem); i++){
		comp = comp*10 + (elem[i] - '0');
	}


	init_CoreHD(argv[1]);
	init_Degree(argv[1]);
	init_Betweenness(argv[1]);
	init_CI(argv[1],rad);

	F = fopen("removedNodes_CoreHD.txt","r"); // metrica 1
	G = fopen("removedNodes_Degree.txt","r"); // metrica 2
	H = fopen("removedNodes_Betweenness.txt","r"); // metrica 3
	I = fopen("removedNodes_CI.txt","r"); // metrica 4
	R = fopen("Kendall_Correlation.csv","a");

	fputs(argv[1],R);
	putc(',',R);
	
	fprintf(stderr, "Calculo metrica CoreHD y Degree\n");

	CalculateKendallMetric(F,G,R,"CoreHD_Degree,",comp);
	putc(',',R);

	fprintf(stderr, "Calculo metrica CoreHD y Betweenness\n");
	
	CalculateKendallMetric(F,H,R,"CoreHD_Betweenness,",comp);
	putc(',',R);

	fprintf(stderr, "Calculo metrica CoreHD y CI\n");
	
	CalculateKendallMetric(F,I,R,"CoreHD_CI,",comp);
	putc(',',R);
		
	fprintf(stderr, "Calculo metrica Degree y Betweenness\n");

	CalculateKendallMetric(G,H,R,"Degree_Betweenness,",comp);
	putc(',',R);

	fprintf(stderr, "Calculo metrica Degree y CI\n");

	CalculateKendallMetric(G,I,R,"Degree_CI,",comp);
	putc(',',R);

	fprintf(stderr, "Calculo metrica Betweenness y CI\n");

	CalculateKendallMetric(H,I,R,"Betweenness_CI,",comp);

	putc('\n',R);

	fclose(R);
	fclose(F);
	fclose(G);
	fclose(H);
	fclose(I);

	return 0;
}
