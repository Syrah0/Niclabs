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


#include "mes.hpp"
#include "omp_config.h"

#include <igraph.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef Mes_t<2> Mes;

struct ConvexCavity {
	void reset() {
		H.clear();
		m1 = INFINITY;
		m2 = INFINITY;
		idxm1 = 0;
		SH = 0;
	}
	void push_back(float h0, float h1) {
		H.push_back(h1);
		SH += h1;
		float h01 = h0 - h1;
		if (h01 <= m1) {
			m2 = m1;
			m1 = h01;
			idxm1 = H.size() - 1;
		} else if (h01 < m2) {
			m2 = h01;
		}

	}
	float cavityval(int i) {
		return SH - H[i];
	}
	float cavitymax(int i) {
		return SH - H[i] + min(0, i == idxm1 ? m2 : m1);
	}
	float fullmax() {
		return SH + min(0, m1);
	}
	float *H;
	float m1, m2, SH;
	int idxm1;
};

struct Buffers {
	Buffers() : Ui(depth + 1), out(depth + 1), C(depth + 1), L0(depth + 1), G1(depth + 1) {}
	void init(int n) {
		//Hin.resize(n, Mes(depth + 1));
		//std::fill(&minh[0], &minh[0] + minh.size(), INFINITY);
		//std::fill(&minh2[0], &minh2[0] + minh2.size(), INFINITY);
		//minh.resize(n, INFINITY);
		//minh2.resize(n, INFINITY);
		for (int t = 0; t <= depth; ++t)
			C[t].reset();
	}
	//Mes *Hin; // tendra tamaño n y cada espacio se inicializa con un mensaje
	//float *minh;
	//float *minh2;
	//Proba Ui;
	//Mes out; //SOLO CONTENDRA DOS PROBAS H[0] Y H[1] (POR TAMAÑO DE n DADO AL TEMPLATE) PARA EL MENSAJE
	// CADA PROBA TENDRA TAMAÑO DEPTH+1
	ConvexCavity *C;
	//float *L0, *G1;
	check_type ck;
};

//Buffers *Mem;

float **Mem = NULL;
float **out = NULL;
float Ui[depth+1];
float L0[depth+1];
float G1[depth+1];

/*struct EdgeProperties {
	EdgeProperties() :
		ij(depth + 1), ji(depth + 1)
	//edge properties
	//messages
	Mes ij, ji;
};
*/
/*
struct VertexProperties  {
	VertexProperties() : H(depth + 1, 0), extH(depth + 1, 0), t(0) {}
	//field
	//Proba H;
	//Proba extH;
	//decisional variable
	//int t;
};*/

/*
typedef adjacency_list<vecS, vecS, undirectedS,
	VertexProperties, EdgeProperties> Graph;
typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
typedef graph_traits<Graph>::out_edge_iterator edge_iterator;
typedef graph_traits<Graph>::edge_iterator graph_edge_iterator;
typedef graph_traits<Graph>::edge_descriptor Edge;
typedef graph_traits<Graph>::vertex_descriptor Vertex;

Graph g;
*/

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
		igraph_neighboors(&graph,&neighbors,u,IGRAPH_OUT);

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
		//g[i].t = times[i]; // arreglar
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

int id_min_element(float *H){
	float min = H[0];
	int t_min = 0;
	int len = sizeof(H)/sizeof(H[0])
	for(int i=1; i < len; i++){
		if(H[i] < min){
			min = H[i];
			t_min = i;
		}
	}
	return t_min;
}

check_type check_v(bool output = false){
	for (int p = 0; p < Mem.size(); ++p)
		Mem[p].ck = check_type();
//#pragma omp parallel for
	for (int j = 0; j < igraph_vcount(&graph); ++j) {
		int sp = 0;
		int tj = igraph_cattribute_VAN(&graph,"t",j);

		igraph_vector_t neighbors, degrees;
		igraph_vector_init(&neighbors,0);
		igraph_neighboors(&graph,&neighbors,vertex,IGRAPH_OUT);

		for (int i = 0; i < igraph_vector_size(&neighbors); i++) {
			int ti = igraph_cattribute_VAN(&graph,"t",igraph_vector_e(&neighbors,i));
			sp += (ti <= tj - 1);
		}
		check_type & ck = Mem[1].ck;
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
		igraph_degree(&graph,&degrees,igraph_vss_all,IGRAPH_OUT,IGRAPH_LOOPS);

		int th = igraph_vector_e(&degrees,j) - 1;
		int good = (th == 0 && tj == 1)
			|| tj == 0
			|| tj == depth
			|| sp >= th;
		if (!good)
			ck.num_bad ++;
	}
	check_type ck;
	for (unsigned p = 0; p < Mem.size(); ++p)
		ck += Mem[p].ck;
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

void free_matrix(void **M){
	int row_Mem = sizeof(M)/sizeof(M[0]);
	for(int i=0; i < row_Mem; i++){
		free(M[i]);
	}
	free(M);
	M = NULL;
}

void init(float **M, int size){
	if(M == NULL){ //init M
		M = (float **)malloc((2*size)*sizeof(float*));

		for(int i = 0; i < (2*size); i++){
			M[i] = (float*)malloc((depth+1)*sizeof(float));
		}
	}
}

float minumum() { // ver cambio a minimo
	float m = inf;
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
	igraph_vector_t result;
	igraph_vector_init(&result, 0);
	igraph_degree(&graph, &result, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS); 
	int n = igraph_vector_e(&result,vertex); // grado del nodo i

	if (n == 0) {
		//g[vertex].H = Proba(depth + 1, INFINITY); // VER
		//g[vertex].H[1] = 0; // g[vertex].H es PROBA -> g[vertex].H[1] = p_[1]
		for(int i=0; i <= depth; i++){
			//igraph_cattribute_VAN_set(&graph,"t",i,ti); // attr t del vertice i
			char name_atr[32];
			sprintf(name_atr,"pdH_%d",i);
			igraph_cattribute_VAN_set(&graph,name_atr,vertex,INFINITY);
		}
		igraph_cattribute_VAN_set(&graph,"pdH_1",vertex,0);
		return 0;
	}
	// siempre se trabaja en base al mismo buffer pq Mem es global y ese se modifica
/*	Buffers & buffers = Mem[0];
	buffers.init(n); 
	VER!!!! ConvexCavity *C = buffers.C;
	LISTO Mes *Hin = buffers.Hin; // tiene n mensajes donde cada uno tiene 2 PROBA de tamaño depth+1 cada uno
	LISTO float *maxh = buffers.maxh;
	LISTO! float *maxh2 = buffers.maxh2;
	LISTO Mes & out = buffers.out; // out tiene 2 PROBA de tamaño depth+1 cada uno
	LISTO!!! Proba & Ui  = buffers.Ui; //ARREGLAR SER USA AL FINAL -> VER SI TRANSFORMAR SOLO A VECTOR
*/

	/*if(Mem == NULL){ //init Mem
		Mem = (float **)malloc((n*2)*sizeof(float*));

		for(int i = 0; i < (2*n); i++){
			Mem[i] = (float*)malloc((depth+1)*sizeof(float));
		}
	}*/
	
	init(Mem, n); // para caso inicial
	init(out, 1); // para caso inicial
	//init_vector(Ui); // para caso inicial

	float Hin[2*n][depth+1];
	float maxh[n];
	float maxh2[n];

	for(int i=0; i<n; i++){
		maxh[i] = INFINITY;
		maxh2[i] = INFINITY;
	}

	Hin = copy_matrix(Mem,Hin);

	free_matrix(Mem); // se libera dado que el siguiente Mem puede tener otras dimensiones

	init(Mem, n); // Se crea Mem con dimensiones correctas

	float summaxh = 0;
	float summaxh2 = 0;

	igraph_vector_t neighbors;
	igraph_vector_init(&neighbors,0);
	igraph_neighboors(&graph,&neighbors,vertex,IGRAPH_OUT);
	int j = 0;

	// igraph_vector_e(&neighbors,k) = vecino k de vertex
	for (int k=0; k < igraph_vector_size(&neighbors); k++, ++j) {
		int index = 2*j; // el otro sera 2*j + 1;
		if(vertex < igraph_vector_e(&neighbors,k)){
			//Hin[j] = g[e].ji; //ARREGLAR 
			//int index = 2*j; // el otro sera 2*j + 1;

			// MES SERA UNA MATRIZ DE 2X(DEPTH+1) TQ MES[0] = H[0] Y MES[1] = H[1]

			igraph_integer_t e;
			igraph_get_eid(&graph,&e,vertex,igraph_vector_e(&neighbors,k),0,0);
			for(int i=0; i<=depth; i++){
				char name_atr[32];
				sprintf(name_atr,"ji0_%d",k);
				Hin[index][i] = igraph_cattribute_EAN(&graph,name_atr,e);
				
				sprintf(name_atr,"ji1_%d",k);
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
				sprintf(name_atr,"ij0_%d",k);
				Hin[index][i] = igraph_cattribute_EAN(&graph,name_atr,e);
				
				sprintf(name_atr,"ij1_%d",k);
				Hin[index+1][i] =  igraph_cattribute_EAN(&graph,name_atr,e);
			}
		}
		//Proba const * h = Hin[j].H; // VA A SER UN VECTOR DE LARGO 2, CONTENDRA AMBOS H DEL MENSAJE Hin
			//AQUIIIII REVISAR!!!!!!!!!!!!
		float *h0 = Hin[index];
		float *h1 = Hin[index+1];

		//float *L0 = buffers.L0; //GLOBALES L0 Y G1!!!!!!!!!!!!!!!!!!!!!!!
		//float *G1 = buffers.G1;
		L0[0] = h0[0]; // h[0] es Hin[j].H[0] (PROBA 0 DEL ARRAY DE 2 PROBAS DEL MENSAJE)
		// Y h[0][0] corresponde al valor almacenado en la p_[0] de PROBA Hin[j].H[0]
		for (int ti = 1; ti <= depth; ++ti){
			//L0[ti] = min(L0[ti - 1],h[0][ti]); // QUIERO EL MINIMO!!
			L0[ti] = min(L0[ti - 1],h0[ti]); // QUIERO EL MINIMO!!
		}
		G1[depth] = INFINITY;
		// REVISAR COMO ACTUA EL ti PQ POR CADA ti 
		for (int ti = depth; ti >= 0; --ti){
			//G1[ti] = min(G1[ti + 1],h[1][ti + 1]); // VER PQ ti + 1 Y NO ti
			G1[ti] = min(G1[ti + 1],h1[ti + 1]); // VER PQ ti + 1 Y NO ti
		}

		for (int ti = 1; ti <= depth; ++ti) {
			float lk = L0[ti - 1]; // Lk
//			float rk = min(h[0][ti],G1[ti]); // VER SI ES MIN O MAX -> SE CONSIDERO MIN
			float rk = min(h0[ti],G1[ti]); // VER SI ES MIN O MAX -> SE CONSIDERO MIN
			C[ti].push_back(rk, lk); // permite calcular M1, M2, k1 VERRRRRRRR
		}

		/* VER ESTAS DEFINICIONES */
//		minh[j] = min(h[0][0], G1[0]); // VER SI ES MIN O MAX --> SE CONSIDERO MIN
		minh[j] = min(h0[0], G1[0]); // VER SI ES MIN O MAX --> SE CONSIDERO MIN
		minh2[j] = L0[depth];
		summinh += minh[j];
		summinh2 += minh2[j];
	}
	//Proba & Hi = g[i].H; // VER COMO AGREGAR ESTE ATRIBUTO AL VERTICE I
	float Hi[depth+1];
	float extH[depth+1];
	for(int k = 0; k <= depth; k++){
		char name_atr[32];
		char name_atr_ext[32];
		sprintf(name_atr,"pdH_%d",k);
		sprintf(name_atr_ext,"pdextH_%d",k);
		Hi[k] = (float)(igraph_cattribute_VAN(&graph,name_atr,i));
		extH[k] = (float)(igraph_cattribute_VAN(&graph,name_atr_ext,i));
	}
	for(int k = 0; k <= depth; k++){
		Hi[k] *= rein;
		Hi[k] += extH[k];
	}

	//Hi *= rein; // multiplica cada elemento de Hi por rein
	//Hi += g[i].extH; // VER COMO AGREGAR ESTE ATRIBUTO AL VERTICE I
	float eps = 0;
	j = 0;
	igraph_vector_destroy(&neighbors);
	
	igraph_vector_init(&neighbors,0);
	igraph_neighboors(&graph,&neighbors,vertex,IGRAPH_OUT);
	for (int k=0; k < igraph_vector_size(&neighbors); k++) {
		int index = 2*j;
		float csumminh = summinh - minh[j]; // VER PARA QUE SIRVE
		float U[2][depth+1];

		U = copy_matrix(out,U);

		//free_matrix(out);
		//Proba * U = out.H; // out es un mensaje => Contine 2 proba (0 y 1) cada una de largo depth+1
		// U contendra las dos proba de out (0 y 1) cada una de largo depth+1

		// OUT ES GLOBAL, U CONTIENE LO QUE CONTIENE OUT Y SI U SE MODIFICA OUT TAMBIEN

		// case ti = 0
		U[0][0] = csumminh - mu; // h0_ij(0)
		U[1][0] = INFINITY; // h1_ij(0)

		// case 0 < ti <= d
		for (int ti = 1; ti < depth; ++ti) {
			U[0][ti] = C[ti].cavityval(j) - ti * rho; // ARREGLAR CAVITY
			U[1][ti] = C[ti].cavitymax(j) - ti * rho; // ARREGLAR CAVITY
		}

		U[0][depth] = summinh2 - minh2[j] - 1;
		U[1][depth] = summinh2 - minh2[j] - 1;

		// A CADA VALOR DE PROBA (0, 1, ..., DEPTH) LE SUMA SU RESPECTIVO DE Hi
		// RECORDAR QUE SE CAMBIO Hi A VECTOR MAS ARRIBA
/*		U[0] += Hi; // h0_ij PROBA 0
		U[1] += Hi; // h1_ij PROBA 1
		out.reduce(); //cambio
*/
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
		if(vertex < igraph_vector_e(&neighbors,k)){
			// AL MENSAJE g[e].ij SE LE ASOCIA EL MENSAJE out (PQ *out Y NO SOLO out???)
			// => g[e].ij TENDRA ASOCIADO A H[O] EL H[0] DE out IDEM PARA H[1]
			igraph_integer_t e;
			igraph_get_eid(&graph,&e,vertex,igraph_vector_e(&neighbors,k),0,0);
			
			//g[e].ij = *out; //ARREGLAR AGREGAR ATTR A LADO
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
			igraph_get_eid(&graph,&e,igraph_vector_e(&neighbors,k),vertex,0,0);

			//g[e].ji = *out; //ARREGLAR AGREGAR ATTR A LADO
			for(int m = 0; m <= depth; m++){
				char name_atr[32];
				sprintf(name_atr,"ji0_%d",k);
				igraph_cattribute_EAN_set(&graph,name_atr,e,out[0][m]);
				
				sprintf(name_atr,"ji1_%d",k);
				igraph_cattribute_EAN_set(&graph,name_atr,e,out[1][m]);
			}
		}
	}

	Ui[0] = summinh - mu;
	for (int ti = 1; ti < depth; ++ti){
		Ui[ti] = C[ti].fullmax() - ti * rho; // verrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
	}
	Ui[depth] =  summinh2 - 1;
	//Hi += Ui;
	for(int i = 0; i <= depth; i++){
		Hi[i] += Ui[i];
	}

	// copio lo de Hin en Mem
	Mem = copy_matrix(Hin,Mem);

	return eps;
}

void converge(){
	int ite = 0;
	float err = 0.0;
	int dec_ite = 0;
	while (err > tolerance && ++ite < maxit && dec_ite < maxdec){
		rein = beta * ite; // tau*gamma
		err = 0;
		int ng = igraph_vcount(&graph); // g grafo, ng numero de vertices 
		for (int i = 0; i < ng; ++i) {
			float diff = update(i); // VER REAL_T
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
	}

	//if (plotting)
	propagate();
}


int main(){
	int depth = 20;
	int maxit = 10000;
	int maxdec = 30;
	float tolerance = 1e-5;
	float beta = 1e-3;
	float rein = 0;
	float noise = 1e-7;
	float mu = 0.1;
	int plotting = false;
	float rho = 1e-5;

	FILE *F;
	char filename[32];
	igraph_t graph;
	igraph_vector_t degrees, nodes, CIvalues;
	double remove = 0.1; // multiplicador de porcentaje
	double rem_nodes = 0.0; // cantidad de nodos removidos
	int total_nodes; // total de nodos del grafo original
	//int T = 35; // limite de la iteracion
	clock_t start, end;
	double time_used;

	G = fopen("MinSum_times.csv","w"); // archivo que guardara los tiempo de ejecucion por iteracion
	H = fopen("MinSum_iter.csv", "w"); // archivo que guardara los R-index (componente mas grande) por iteracion

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

	for(int i=0; i < igraph_vcount(&graph); i++){
		igraph_vector_push_back(&attr_t_vertex,0);
		igraph_vector_push_back(&attr_depth_vertex,depth+1);
	}
	igraph_cattribute_VAN_setv(&graph, "t", &attr_t_vertex);
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

	// check_v(true); // VER

	igraph_vector_t rem_vertex;
	igraph_vs_t delete_vertex;
	igraph_vector_init(&rem_vertex);
	for(int j = 0; j < igraph_vcount(&graph); j++){
		/*if(g[j].t == 0){ //hay que removerlo VERRR
			igraph_vector_push_back(&rem_vertex, j);
		}*/
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
