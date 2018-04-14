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
	return res/div;
}

double CalculateKendallMetric(FILE *F, FILE *G, int comp){
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

	//fprintf(stderr, "%lf\n", kendall);
	return kendall;
}

char** initKendall(char *argv[], int pos){

	FILE *F, *G, *H, *I, *R;
	int rad = atoi(argv[1]);
	int comp = atoi(argv[3]);
	char filename[300];
	char value[32];
	char **ret;
	ret = malloc(12 * sizeof(char *));
	for(int i = 0; i < 12; i++){
		ret[i] = malloc(50);
	}
	double res;

	sprintf(filename, "archivos/CoreHD/removedNodes_CoreHD_graph%d.txt",pos);

	init_CoreHD(argv[0],comp, pos);

	F = fopen(filename,"r"); // metrica 1
	
	sprintf(filename, "archivos/Degree/removedNodes_Degree_graph%d.txt",pos);

	init_Degree(argv[0],comp, pos);
	
	G = fopen(filename,"r"); // metrica 2

	sprintf(filename, "archivos/Betweenness/removedNodes_Betweenness_graph%d.txt",pos);
	
	init_Betweenness(argv[0],comp, pos);
	
	H = fopen(filename,"r"); // metrica 3
	
	sprintf(filename, "archivos/CI/removedNodes_CI_graph%d.txt",pos);
	
	init_CI(argv[0],rad,argv[2],comp, pos);
	
	I = fopen(filename,"r"); // metrica 4
	
	fprintf(stderr, "Corriendo KENDALL\n");
	fprintf(stderr, "\n");
	
	fprintf(stderr, "Calculo metrica CoreHD y Degree\n");

	res = CalculateKendallMetric(F,G,comp);

	sprintf(ret[1],"%f",res);
	ret[0] = "CoreHD_Degree,";

	fprintf(stderr, "Calculo metrica CoreHD y Betweenness\n");
	
	res = CalculateKendallMetric(F,H,comp);

	sprintf(ret[3],"%f",res);
	ret[2] = "CoreHD_Betweenness,";

	fprintf(stderr, "Calculo metrica CoreHD y CI\n");
	
	res = CalculateKendallMetric(F,I,comp);
	
	sprintf(ret[5],"%f",res);
	ret[4] = "CoreHD_CI,";
		
	fprintf(stderr, "Calculo metrica Degree y Betweenness\n");

	res = CalculateKendallMetric(G,H,comp);
	
	sprintf(ret[7],"%f",res);
	ret[6] = "Degree_Betweenness,";

	fprintf(stderr, "Calculo metrica Degree y CI\n");

	res = CalculateKendallMetric(G,I,comp);
	
	sprintf(ret[9],"%f",res);
	ret[8] = "Degree_CI,";
	fprintf(stderr, "Calculo metrica Betweenness y CI\n");

	res = CalculateKendallMetric(H,I,comp);
	
	sprintf(ret[11],"%f",res);
	ret[10] = "Betweenness_CI,";

	fclose(F);
	fclose(G);
	fclose(H);
	fclose(I);

	return ret;
}
