/*
 * Decycler, a reinforced Max-Sum algorithm to solve the decycling problem
 * Copyright (C) 2016 Alfredo Braunstein
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 2 of the License.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/


#include <igraph.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


struct ConvexCavity {
	igraph_vector_t H;
	float m1, m2, SH;
	int idxm1;

	//void (*reset)(struct ConvexCavity *);
	/*void (*push_back)(struct ConvexCavity *, float h0, float h1);
	float (*cavityval)(struct ConvexCavity *, int i);
	float (*cavitymax)(struct ConvexCavity *, int i);
	float (*fullmax)(struct ConvexCavity *);*/
};

typedef struct{
	int num_bad;
	int num_seeds;
	int num_on;
	float tot_cost;
	float tot_energy;
} check_type;


#define depth 20
#define maxit 10000
#define maxdec 30
#define tolerance 1e-5
#define beta 1e-3
#define noise 1e-7
#define mu 0.1
#define rho 1e-5
int rein = 0;
float **Mem = NULL;
float **out = NULL;
float Ui[depth+1];
float L0[depth+1];
float G1[depth+1];
igraph_t graph;

void propagate(){
	int theta[igraph_vcount(&graph)];
	int times[igraph_vcount(&graph)];
	
	igraph_vector_t degrees;
	igraph_vector_init(&degrees,0);
	igraph_degree(&graph, &degrees, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS); 

	igraph_dqueue_t q;
	igraph_dqueue_init(&q,0);
	
	for(int v = 0; v < igraph_vcount(&graph); v++) {
		int val_t = (int)(igraph_cattribute_VAN(&graph,"t",v));
		if (val_t == 0) {
			theta[v] = 0;
			times[v] = 0;
			igraph_dqueue_push(&q, v);
		} else {
			theta[v] = igraph_vector_e(&degrees,v) - 1;
			times[v] = depth;
			if (theta[v]) {
				times[v] = 1;
				igraph_dqueue_push(&q, v);
			}
		}
	}
	while (!igraph_dqueue_empty(&q)) {
		int u = (int)igraph_dqueue_pop(&q);
		//edge_iterator eit, end;
		igraph_vector_t neighbors;
		igraph_vector_init(&neighbors,0);
		igraph_neighbors(&graph,&neighbors,u,IGRAPH_OUT);

		for (int i = 0; i < igraph_vector_size(&neighbors); i++) {
			int v = igraph_vector_e(&neighbors,i);
			if (--theta[v] == 0) {
				times[v] = times[u] + 1;
				igraph_dqueue_push(&q, v);
			}
		}
	}
	int seeds = 0;
	int num_on = 0;
	for (int i = 0; i < igraph_vcount(&graph); i++) {
		num_on += (times[i] < depth);
		seeds += (times[i] == 0);
		int new_val = times[i];
		igraph_cattribute_VAN_set(&graph,"t",i,new_val);
	}
}

int min(int a, int b){
	if(a>b){
		return b;
	}
	else{
		return a;
	}
}

int max(int a, int b){
	if(a<b){
		return b;
	}
	else{
		return a;
	}
}

struct ConvexCavity reset(struct ConvexCavity convex) {
	//fprintf(stderr, "AQUIII\n");
	igraph_vector_init(&(convex.H),0);
	convex.m1 = INFINITY;
	convex.m2 = INFINITY;
	convex.idxm1 = 0;
	convex.SH = 0;
	return convex;
}

struct ConvexCavity push_back(struct ConvexCavity convex,float h0, float h1) {
	igraph_vector_push_back(&(convex.H),h1);
	convex.SH += h1;
	float h01 = h0 - h1;
	if (h01 <= (convex.m1)) {
		convex.m2 = convex.m1;
		convex.m1 = h01;
		convex.idxm1 = ((int)igraph_vector_size(&(convex.H))) - 1;
	} 
	else if (h01 < (convex.m2)) {
		convex.m2 = h01;
	}
	return convex;
}

float cavityval(struct ConvexCavity convex,int i) {
	return (convex.SH) - ((float)igraph_vector_e(&(convex.H),i));
}

float cavitymax(struct ConvexCavity convex,int i) {
	return (convex.SH) - ((float)igraph_vector_e(&(convex.H),i)) + min(0, i == (convex.idxm1) ? (convex.m2) : (convex.m1));
}

float fullmax(struct ConvexCavity convex) {
	return (convex.SH) + min(0, (convex.m1));
}

int id_min_element(float *H){
	float min = H[0];
	int t_min = 0;
	int len = sizeof(H)/sizeof(H[0]);
	for(int i=1; i < len; i++){
		if(H[i] < min){
			min = H[i];
			t_min = i;
		}
	}
	return t_min;
}

check_type check_v() {
	check_type ck = {0,0,0,0,0};
	for (int j = 0; j < igraph_vcount(&graph); ++j) {
		igraph_vector_t neighbors, degrees;
		int sp = 0;
		int tj = (int) igraph_cattribute_VAN(&graph,"t",j);

		igraph_vector_init(&neighbors,0);
		igraph_neighbors(&graph,&neighbors,j,IGRAPH_OUT);

		// igraph_vector_e(&neighbors,k) = vecino k de vertex
		for (int k=0; k < igraph_vector_size(&neighbors); k++) {
			int ti = (int) igraph_cattribute_VAN(&graph,"t",(igraph_vector_e(&neighbors,k)));
			sp += (ti <= tj - 1);
		}
		if (tj < depth) {
			ck.num_on++;
			ck.tot_energy--;
		}
		if (tj == 0) {
			ck.num_seeds++;
			ck.tot_cost++;
			ck.tot_energy += mu;
		}
		igraph_vector_init(&degrees,0);
		igraph_degree(&graph, &degrees, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
		
		int th = igraph_vector_e(&degrees,j) - 1;
		int good = (th == 0 && tj == 1)
			|| tj == 0
			|| tj == depth
			|| sp >= th;
		if (!good)
			ck.num_bad ++;

		igraph_vector_destroy(&neighbors);
		igraph_vector_destroy(&degrees);
	}
	return ck;
}

float **copy_matrix(float **M, float **Hin){


	int row_Mem = sizeof(M)/sizeof(M[0]);
	int row_Hin = sizeof(Hin)/sizeof(Hin[0]);

	if(row_Mem <= row_Hin){
		for(int i=0; i < row_Mem; i++){
			for(int j=0; j <= depth; j++){
				Hin[i][j] = M[i][j];
			}
		}
		if(row_Mem != row_Hin){ //rellenar el resto de Hin con 0's
			for(int i=row_Mem; i < row_Hin; i++){
				for(int j=0; j <= depth; j++){
					Hin[i][j] = 0;
				}
			}
		}
	}
	else{
		for(int i=0; i < row_Hin; i++){
			for(int j=0; j <= depth; j++){
				Hin[i][j] = M[i][j];
			}
		}
	}

	return Hin;
}

float** free_matrix(float **M){
	int row_Mem = sizeof(M)/sizeof(M[0]);
	for(int i=0; i < row_Mem; i++){
		free(M[i]);
	}
	free(M);
	M = NULL;
	return M;
}

float** init(float **M, int size){
	if(M == NULL){ //init M
		M = (float **)malloc((2*size)*sizeof(float*));

		for(int i = 0; i < (2*size); i++){
			M[i] = (float*)malloc((depth+1)*sizeof(float));
		}

		for(int i=0; i < (2*size); i++){
			for(int j=0; j <= depth; j++){
				M[i][j] = 0;
			}
		}
	}
	return M;
}

float minimum() { // ver cambio a minimo
	float m = INFINITY;
	//int const depth = size();
	for (int i = 0; i < 2; i++)
		for (int j = 0; j <= depth; j++)
			m = min(m, out[i][j]);
	return m;
}

void reduce(){
	float m = minimum();
	//int const depth = size();
	for (int i = 0; i < 2; i++)
		for (int j = 0; j <= depth; j++)
			out[i][j] -= m; // revisar si es resta
}

float l8dist_aux(float *old, float *out){
	float n = 0;
	for(int i = depth+1; i >= 0; i--){
		n = max(n, fabs(old[i] - out[i]));
	}
	return n;
}

float l8dist(float *old0, float *old1, float **out){ 
	float l8 = 0;
	l8 = max(l8,l8dist_aux(old0,out[0]));
	l8 = max(l8,l8dist_aux(old1,out[1]));
	return l8;
}

float update(int vertex){ // i = numero que representa al vertice.{
	igraph_vector_t degrees;
	igraph_vector_init(&degrees, 0);
	igraph_degree(&graph, &degrees, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS); 
	int n = igraph_vector_e(&degrees,vertex); // grado del nodo i


	if (n == 0) {
		for(int i=0; i <= depth; i++){
			char name_atr[32];
			sprintf(name_atr,"pdH_%d",i);
			igraph_cattribute_VAN_set(&graph,name_atr,vertex,INFINITY);
		}
		igraph_cattribute_VAN_set(&graph,"pdH_1",vertex,0);
		return 0;
	}
	// siempre se trabaja en base al mismo buffer pq Mem es global y ese se modifica

	struct ConvexCavity C[depth+1]; 
	for(int i = 0; i<=depth; i++){
		C[i] = reset(C[i]);
	}
	
	Mem = init(Mem, n); // para caso inicial
	out = init(out, 1); // para caso inicial
	
	float **Hin = NULL;
	Hin = init(Hin, n);
	
	float minh[n];
	float minh2[n];

	for(int i=0; i<n; i++){
		minh[i] = INFINITY;
		minh2[i] = INFINITY;
	}

	Hin = copy_matrix(Mem,Hin);

	Mem = free_matrix(Mem); // se libera dado que el siguiente Mem puede tener otras dimensiones

	Mem = init(Mem, n); // Se crea Mem con dimensiones correctas

	float summinh = 0;
	float summinh2 = 0;

	igraph_vector_t neighbors;
	igraph_vector_init(&neighbors,0);
	igraph_neighbors(&graph,&neighbors,vertex,IGRAPH_OUT);
	int j = 0;

	// igraph_vector_e(&neighbors,k) = vecino k de vertex
	for (int k=0; k < igraph_vector_size(&neighbors); k++, ++j) {
		int index = 2*j; // el otro sera 2*j + 1;
		if(vertex < igraph_vector_e(&neighbors,k)){

			// MES SERA UNA MATRIZ DE 2X(DEPTH+1) TQ MES[0] = H[0] Y MES[1] = H[1]

			igraph_integer_t e;
			igraph_get_eid(&graph,&e,vertex,igraph_vector_e(&neighbors,k),0,0);
			for(int i=0; i<=depth; i++){
				char name_atr[32];
				sprintf(name_atr,"ji0_%d",i);
				Hin[index][i] = igraph_cattribute_EAN(&graph,name_atr,e);
				
				sprintf(name_atr,"ji1_%d",i);
				Hin[index+1][i] =  igraph_cattribute_EAN(&graph,name_atr,e);
			}
		}
		else{
			// Hin[j] = Corresponde al j-esimo mensaje del vector Hin
			// este mensaje contiene 2 Proba de largo depth+1 cada una
			// => le asocia a g[e].ij las dos Proba de Hin[j]
			//Hin[j] = g[e].ij; //ARREGLAR	

			igraph_integer_t e;
			igraph_get_eid(&graph,&e,igraph_vector_e(&neighbors,k),vertex,0,0);
			// MES SERA UNA MATRIZ DE 2X(DEPTH+1) TQ MES[0] = H[0] Y MES[1] = H[1]
			for(int i=0; i<=depth; i++){
				char name_atr[32];
				sprintf(name_atr,"ij0_%d",i);
				Hin[index][i] = igraph_cattribute_EAN(&graph,name_atr,e);
				
				sprintf(name_atr,"ij1_%d",i);
				Hin[index+1][i] =  igraph_cattribute_EAN(&graph,name_atr,e);
			}
		}
		//Proba const * h = Hin[j].H; // VA A SER UN VECTOR DE LARGO 2, CONTENDRA AMBOS H DEL MENSAJE Hin
		float *h0 = Hin[index];
		float *h1 = Hin[index+1];

		L0[0] = h0[0]; // h[0] es Hin[j].H[0] (PROBA 0 DEL ARRAY DE 2 PROBAS DEL MENSAJE)
		// Y h[0][0] corresponde al valor almacenado en la p_[0] de PROBA Hin[j].H[0]
		for (int ti = 1; ti <= depth; ++ti){
			//L0[ti] = min(L0[ti - 1],h[0][ti]); // QUIERO EL MINIMO!!
			L0[ti] = min(L0[ti - 1],h0[ti]); // QUIERO EL MINIMO!!
		}
		G1[depth] = INFINITY;
		// REVISAR COMO ACTUA EL ti PQ POR CADA ti 
		for (int ti = depth; ti >= 0; --ti){
			G1[ti] = min(G1[ti + 1],h1[ti + 1]); // VER PQ ti + 1 Y NO ti
		}

		for (int ti = 1; ti <= depth; ++ti) {
			float lk = L0[ti - 1]; // Lk
			float rk = min(h0[ti],G1[ti]); // VER SI ES MIN O MAX -> SE CONSIDERO MIN
			C[ti] = push_back(C[ti],rk, lk); // permite calcular M1, M2, k1 VERRRRRRRR
		}

		/* VER ESTAS DEFINICIONES */
		minh[j] = min(h0[0], G1[0]); // VER SI ES MIN O MAX --> SE CONSIDERO MIN
		minh2[j] = L0[depth];
		summinh += minh[j];
		summinh2 += minh2[j];
	}
	float Hi[depth+1];
	float extH[depth+1];
	for(int k = 0; k <= depth; k++){
		char name_atr[32];
		char name_atr_ext[32];
		sprintf(name_atr,"pdH_%d",k);
		sprintf(name_atr_ext,"pdextH_%d",k);
		Hi[k] = (float)(igraph_cattribute_VAN(&graph,name_atr,vertex));
		extH[k] = (float)(igraph_cattribute_VAN(&graph,name_atr_ext,vertex));
	}
	for(int k = 0; k <= depth; k++){
		Hi[k] *= rein;
		Hi[k] += extH[k];
	}

	float eps = 0;
	j = 0;
	igraph_vector_destroy(&neighbors);
	
	igraph_vector_init(&neighbors,0);
	igraph_neighbors(&graph,&neighbors,vertex,IGRAPH_OUT);
	int size = igraph_vector_size(&neighbors);
	for (int k=0; k < size; k++) {
		//fprintf(stderr, "Aqui\n");
		int index = 2*j;
		int neigh = (int)igraph_vector_e(&neighbors,k);
		float csumminh = summinh - minh[j]; // VER PARA QUE SIRVE
		float **U = NULL;
		U = init(U,1);
		U = copy_matrix(out,U);

		//Proba * U = out.H; // out es un mensaje => Contine 2 proba (0 y 1) cada una de largo depth+1
		// U contendra las dos proba de out (0 y 1) cada una de largo depth+1

		// OUT ES GLOBAL, U CONTIENE LO QUE CONTIENE OUT Y SI U SE MODIFICA OUT TAMBIEN

		// case ti = 0
		U[0][0] = csumminh - mu; // h0_ij(0)
		U[1][0] = INFINITY; // h1_ij(0)

		// case 0 < ti <= d
		for (int ti = 1; ti < depth; ++ti) {
			U[0][ti] = cavityval(C[ti],j) - ti * rho; // ARREGLAR CAVITY
			U[1][ti] = cavitymax(C[ti],j) - ti * rho; // ARREGLAR CAVITY
		}

		U[0][depth] = summinh2 - minh2[j] - 1;
		U[1][depth] = summinh2 - minh2[j] - 1;

		// A CADA VALOR DE PROBA (0, 1, ..., DEPTH) LE SUMA SU RESPECTIVO DE Hi
		// RECORDAR QUE SE CAMBIO Hi A VECTOR MAS ARRIBA
		for(int m = 0; m <= depth; m++){
			U[0][m] += Hi[m];
			U[1][m] += Hi[m];
		}

		out = copy_matrix(U,out);
		reduce();

		//Mes & old = Hin[j]; // old SERAN 2 PROBA (0 Y 1) -> old CONTENDRA LAS 2 PROBA ASOCIADAS AL MENSAJE j-esimo DE Hin
		float *old0 =  Hin[index];
		float *old1  = Hin[index+1];

		eps = l8dist(old0, old1, out); //cambio
		// set mes
		// VER QUE GUARDA EL MENSAJE QUE SE PASA PARA PODER DECODIFICARLO
		if(vertex < neigh){
			// AL MENSAJE g[e].ij SE LE ASOCIA EL MENSAJE out (PQ *out Y NO SOLO out???)
			// => g[e].ij TENDRA ASOCIADO A H[O] EL H[0] DE out IDEM PARA H[1]
			igraph_integer_t e;
			igraph_get_eid(&graph,&e,vertex,neigh,0,0);
			
			for(int m = 0; m <= depth; m++){
				char name_atr[32];
				sprintf(name_atr,"ij0_%d",k);
				igraph_cattribute_EAN_set(&graph,name_atr,e,out[0][m]);
				sprintf(name_atr,"ij1_%d",k);
				igraph_cattribute_EAN_set(&graph,name_atr,e,out[1][m]);
			}
		}
		else{
			igraph_integer_t e;
			igraph_get_eid(&graph,&e,neigh,vertex,0,0);
			for(int m = 0; m <= depth; m++){
				char name_atr[32];
				sprintf(name_atr,"ji0_%d",k);
				igraph_cattribute_EAN_set(&graph,name_atr,e,out[0][m]);
				
				sprintf(name_atr,"ji1_%d",k);
				igraph_cattribute_EAN_set(&graph,name_atr,e,out[1][m]);
			}
		}
		U = free_matrix(U);
	}
	igraph_vector_destroy(&degrees);

	igraph_vector_destroy(&neighbors);
	//fprintf(stderr, "HOLAA\n");
	Ui[0] = summinh - mu;
	for (int ti = 1; ti < depth; ++ti){
		Ui[ti] = fullmax(C[ti]) - ti * rho; // verrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
	}
	Ui[depth] =  summinh2 - 1;
	for(int i = 0; i <= depth; i++){
		Hi[i] += Ui[i];
	}

	// copio lo de Hin en Mem
	Mem = copy_matrix(Hin,Mem);
	Hin = free_matrix(Hin);

	return eps;
}

void converge(){
	int ite = 0;
	float err = 0.0;
	int dec_ite = 0;
	do{
		rein = beta * ite; // tau*gamma
		err = 0;
		int ng = igraph_vcount(&graph); // g grafo, ng numero de vertices 
		for (int i = 0; i < ng; ++i) {
			float diff = update(i); // VER REAL_T
			//fprintf(stderr, "AQUI2\n" );
			if (diff > err){
				err = diff;
			}
		}
		++dec_ite;
		int numon2 = 0;
		for (int i = 0; i < ng; ++i) {
			// g[i].H es proba -> g[i].H[k] poner en Hi[k] (depth sera el largo)
			//Proba & Hi = g[i].H; // VER COMO OBTENER EL ATTR H DEL VERTICE I
			float Hi[depth+1];
			for(int k = 0; k <= depth; k++){
				char name_atr[32];
				sprintf(name_atr,"pdH_%d",k);
				Hi[k] = (float)(igraph_cattribute_VAN(&graph,name_atr,i));
			}
			int ti = id_min_element(Hi);
			float Hmax = Hi[ti];
			int len = sizeof(Hi)/sizeof(Hi[0]);
			for (int t = 0; t < len; ++t) {
				Hi[t] -= Hmax;;
			}
			numon2 += (ti < depth);

			/*if (ti != g[i].t){
				dec_ite = 0;
			}*/
			int val_t = (int)(igraph_cattribute_VAN(&graph,"t",i));
			if (ti != val_t){
				dec_ite = 0;
			}
			//g[i].t = ti; // attr t del vertice i
			igraph_cattribute_VAN_set(&graph,"t",i,ti); // attr t del vertice i
		}
		check_type ck = check_v();
		if (ck.num_bad){
			dec_ite = 0;
		}
		fprintf(stderr, "IT: %d/%d | DEC: %d/%d | ERR: %lf/%lf | ON: %d | S: %d | ON2: %d | BAD: %d | COST: %lf | EN: %lf\n",ite,maxit,dec_ite,maxdec,err,tolerance,ck.num_on,ck.num_seeds,numon2,ck.num_bad,ck.tot_cost,ck.tot_energy);
	}while (err > tolerance && ++ite < maxit && dec_ite < maxdec);
	//if (plotting)
	propagate();
}


int main(){

	FILE *F, *G, *H;
	char filename[32];
	igraph_vector_t degrees, nodes, CIvalues;
	double remove = 0.1; // multiplicador de porcentaje
	double rem_nodes = 0.0; // cantidad de nodos removidos
	int total_nodes; // total de nodos del grafo original
	//int T = 35; // limite de la iteracion
	clock_t start, end;
	double time_used;

	G = fopen("MinSum_times.csv","w"); // archivo que guardara los tiempo de ejecucion por iteracion
	H = fopen("MinSum_iter.csv", "w"); // archivo que guardara los R-index (componente mas grande) por iteracion

	/* turn on attribute handling */
	igraph_i_set_attribute_table(&igraph_cattribute_table);

	/* Se lee el archivo que contiene las conexiones de los nodos */
	F = fopen("red3.edges","r");
	igraph_read_graph_edgelist(&graph,F,0,0); // crea el grafo a partir del archivo con las conexiones
	fclose(F);
	// agregar atributo a todos los vertices del grafo
	// atribute t permitira verificar si debe ser eliminado o no
	igraph_vector_t attr_t_vertex;
	igraph_vector_t attr_depth_vertex;
	igraph_vector_init(&attr_t_vertex, 0);
	igraph_vector_init(&attr_depth_vertex, 0);

	//igraph_vector_fill(&attr_t_vertex,0);
	for(int i=0; i < igraph_vcount(&graph); i++){
		igraph_vector_push_back(&attr_t_vertex,0);
		igraph_vector_push_back(&attr_depth_vertex,depth+1);
	}
	igraph_cattribute_VAN_setv(&graph, "t", &attr_t_vertex);
	//SETVANV(&g, "t", &attr_t_vertex);
	igraph_cattribute_VAN_setv(&graph, "depth", &attr_depth_vertex);

	// set atributos H y extH -> son "depth" atributos
	for(int i=0; i<=depth; i++){
		igraph_vector_t attr_pd_vertex;
		igraph_vector_t attr_pdextH_vertex;
		char name_atr[32];
		char name_atr_ext[32];

		igraph_vector_init(&attr_pd_vertex,0);
		igraph_vector_init(&attr_pdextH_vertex,0);
		
		sprintf(name_atr,"pdH_%d",i);
		sprintf(name_atr_ext,"pdextH_%d",i);

		for(int j=0; j < igraph_vcount(&graph); j++){
			igraph_vector_push_back(&attr_pd_vertex,0);
			igraph_vector_push_back(&attr_pdextH_vertex,0);
		}

		igraph_cattribute_VAN_setv(&graph,name_atr,&attr_pd_vertex);
		igraph_cattribute_VAN_setv(&graph,name_atr_ext,&attr_pdextH_vertex);
		
		igraph_vector_destroy(&attr_pdextH_vertex);
		igraph_vector_destroy(&attr_pd_vertex);
	}

	for(int i=0; i<=depth; i++){
		igraph_vector_t attr_ij0_edge;
		igraph_vector_t attr_ij1_edge;
		igraph_vector_t attr_ji0_edge;
		igraph_vector_t attr_ji1_edge;
		char name_atr_ij0[32];
		char name_atr_ji0[32];
		char name_atr_ij1[32];
		char name_atr_ji1[32];

		igraph_vector_init(&attr_ij0_edge,0);
		igraph_vector_init(&attr_ji0_edge,0);
		igraph_vector_init(&attr_ij1_edge,0);
		igraph_vector_init(&attr_ji1_edge,0);
		
		sprintf(name_atr_ij0,"ij0_%d",i);
		sprintf(name_atr_ji0,"ji0_%d",i);
		sprintf(name_atr_ij1,"ij1_%d",i);
		sprintf(name_atr_ji1,"ji1_%d",i);


		for(int j=0; j < igraph_ecount(&graph); j++){
			igraph_vector_push_back(&attr_ji0_edge,0);
			igraph_vector_push_back(&attr_ij0_edge,0);
			igraph_vector_push_back(&attr_ji1_edge,0);
			igraph_vector_push_back(&attr_ij1_edge,0);
		}

		igraph_cattribute_EAN_setv(&graph,name_atr_ij0,&attr_ij0_edge);
		igraph_cattribute_EAN_setv(&graph,name_atr_ji0,&attr_ji0_edge);
		igraph_cattribute_EAN_setv(&graph,name_atr_ij1,&attr_ij1_edge);
		igraph_cattribute_EAN_setv(&graph,name_atr_ji1,&attr_ji1_edge);
		
		igraph_vector_destroy(&attr_ij1_edge);
		igraph_vector_destroy(&attr_ij0_edge);
		igraph_vector_destroy(&attr_ji1_edge);
		igraph_vector_destroy(&attr_ji0_edge);
	}

	//plotting = true; // VER


	converge(); 

	//fprintf(stderr, "AQUIIIIII\n");
	// check_v(true); // VER

	igraph_vector_t rem_vertex;
	igraph_vs_t delete_vertex;
	igraph_vector_init(&rem_vertex,0);
	for(int j = 0; j < igraph_vcount(&graph); j++){
		int val_t = (int)(igraph_cattribute_VAN(&graph,"t",j));
		if(val_t == 0){ //hay que removerlo VERRR
			igraph_vector_push_back(&rem_vertex, j);
		}
	}

	// SETEAR CARACTERISTICA DE PROBA -> 1 P_[D] POR CADA P_[D] QUE DEBERIA TENER PROBA -> HACER CON FOR ESTO POR CADA NODO
	// SETEAR DEPTH COMO ATTR NUMERICO PROPIO DE CADA VERTICE -> DEPTH NUNCA CAMBIA VER SI SE ALTERA EN PROBA O MES

	igraph_vs_vector(&delete_vertex,&rem_vertex);
	igraph_delete_vertices(&graph,delete_vertex);

	igraph_vector_destroy(&rem_vertex);
	igraph_vs_destroy(&delete_vertex);
	igraph_destroy(&graph);

	//COMO VER EL 10% - 20% ... DE REMOSION

	return 0;
}
