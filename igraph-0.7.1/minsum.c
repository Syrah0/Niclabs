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
		real_t const h01 = h0 - h1;
		if (h01 <= m1) {
			m2 = m1;
			m1 = h01;
			idxm1 = H.size() - 1;
		} else if (h01 < m2) {
			m2 = h01;
		}

	}
	float cavityval(int i) const {
		return SH - H[i];
	}
	float cavitymax(int i) const {
		return SH - H[i] + min(0, i == idxm1 ? m2 : m1);
	}
	float fullmax() const {
		return SH + min(0, m1);
	}
	float *H;
	float m1, m2, SH;
	int idxm1;
};

struct Buffers {
	Buffers() : Ui(depth + 1), out(depth + 1), C(depth + 1), L0(depth + 1), G1(depth + 1) {}
	void init(int n) {
		Hin.resize(n, Mes(depth + 1));
		std::fill(&minh[0], &minh[0] + minh.size(), INFINITY);
		std::fill(&minh2[0], &minh2[0] + minh2.size(), INFINITY);
		minh.resize(n, INFINITY);
		minh2.resize(n, INFINITY);
		for (int t = 0; t <= depth; ++t)
			C[t].reset();
	}
	Mes *Hin;
	float *minh;
	float *minh2;
	Proba Ui;
	Mes out;
	ConvexCavity *C;
	float *L0, *G1;
};

Buffers *Mem;

struct EdgeProperties {
	EdgeProperties() :
		ij(depth + 1), ji(depth + 1)
	//edge properties
	//messages
	Mes ij, ji;
};

struct VertexProperties  {
	VertexProperties() : H(depth + 1, 0), extH(depth + 1, 0), t(0) {}
	//field
	Proba H;
	Proba extH;
	//decisional variable
	//int t;
};

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

int id_min_element(Proba H){
	float min = H[0];
	int t_min = 0;
	for(int i=1; i < (int)(H.size()); i++){
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

float update(int vertex){ // i = numero que representa al vertice.{
	igraph_vector_t result;
	igraph_vector_init(&result, 0);
	igraph_degree(&graph, &result, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS); 
	int n = igraph_vector_e(&result,vertex); // grado del nodo i

	if (n == 0) {
		g[vertex].H = Proba(depth + 1, INFINITY); // VER
		g[vertex].H[1] = 0;
		return 0;
	}

	Buffers & buffers = Mem[1];
	buffers.init(n);
	ConvexCavity *C = buffers.C;
	Mes *Hin = buffers.Hin;
	float *maxh = buffers.maxh;
	float *maxh2 = buffers.maxh2;
	Mes & out = buffers.out;
	Proba & Ui  = buffers.Ui;

	float summaxh = 0;
	float summaxh2 = 0;

	igraph_vector_t neighbors;
	igraph_vector_init(&neighbors,0);
	igraph_neighboors(&graph,&neighbors,vertex,IGRAPH_OUT);
	int j = 0;

	// igraph_vector_e(&neighbors,k) = vecino k de vertex
	for (int k=0; k < igraph_vector_size(&neighbors); k++) {
		if(vertex < igraph_vector_e(&neighbors,k)){
			Hin[j] = g[e].ji; //ARREGLAR 
		}
		else{
			Hin[j] = g[e].ij; //ARREGLAR	
		}
		Proba const * h = Hin[j].H;
		int *L0 = buffers.L0;
		int *G1 = buffers.G1;
		L0[0] = h[0][0];
		for (int ti = 1; ti <= depth; ++ti){
			L0[ti] = min(L0[ti - 1],h[0][ti]); // QUIERO EL MINIMO!!
		}
		G1[depth] = INFINITY;
		// REVISAR COMO ACTUA EL ti PQ POR CADA ti 
		for (int ti = depth; ti >= 0; --ti){
			G1[ti] = min(G1[ti + 1],h[1][ti + 1]); // VER PQ ti + 1 Y NO ti
		}

		for (int ti = 1; ti <= depth; ++ti) {
			float lk = L0[ti - 1]; // Lk
			float rk = min(h[0][ti],G1[ti]); // VER SI ES MIN O MAX -> SE CONSIDERO MIN
			C[ti].push_back(rk, lk); // permite calcular M1, M2, k1 VERRRRRRRR
		}

		/* VER ESTAS DEFINICIONES */
		minh[j] = min(h[0][0], G1[0]); // VER SI ES MIN O MAX --> SE CONSIDERO MIN
		minh2[j] = L0[depth];
		summinh += minh[j];
		summinh2 += minh2[j];
	}
	Proba & Hi = g[i].H; // VER COMO AGREGAR ESTE ATRIBUTO AL VERTICE I
	Hi *= rein; // multiplica cada elemento de Hi por rein
	Hi += g[i].extH; // VER COMO AGREGAR ESTE ATRIBUTO AL VERTICE I
	float eps = 0;
	j = 0;
	igraph_vector_destroy(&neighbors);
	
	igraph_vector_init(&neighbors,0);
	igraph_neighboors(&graph,&neighbors,vertex,IGRAPH_OUT);
	for (int k=0; k < igraph_vector_size(&neighbors); k++) {
		float csumminh = summinh - minh[j]; // VER PARA QUE SIRVE
		Proba * U = out.H;

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

		U[0] += Hi;
		U[1] += Hi;
		out.reduce(); //cambio

		Mes & old = Hin[j];
		eps = l8dist(old, out); //cambio
		// set mes
		if(vertex < igraph_vector_e(&neighbors,k)){
			g[e].ij = *out; //ARREGLAR AGREGAR ATTR A LADO
		}
		else{
			g[e].ji = *out; //ARREGLAR AGREGAR ATTR A LADO
		}
	}

	Ui[0] = summinh - mu;
	for (int ti = 1; ti < depth; ++ti){
		Ui[ti] = C[ti].fullmax() - ti * rho; // verrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
	}
	Ui[depth] =  summinh2 - 1;
	Hi += Ui;
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
///////////AQUIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
		++dec_ite;
		int numon2 = 0;
		for (int i = 0; i < ng; ++i) {
			Proba & Hi = g[i].H; // VER COMO OBTENER EL ATTR H DEL VERTICE I
			int ti = id_min_element(Hi);
			double Hmax = Hi[ti];
			for (int t = 0; t < int(Hi.size()); ++t) {
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
	igraph_vector_init(&attr_t_vertex, 0);

	for(int i=0; i < igraph_vcount(&graph); i++){
		igraph_vector_push_back(&attr_t_vertex,0);
	}
	igraph_cattribute_VAN_setv(&graph, "t", &attr_t_vertex);

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
	// SETEAR DEPTH COMO ATTR NUMERICO PROPIO DE CADA VERTICE 

	igraph_vs_vector(&delete_vertex,&rem_vertex);
	igraph_delete_vertices(&graph,delete_vertex);

	igraph_vector_destroy(&rem_vertex);
	igraph_vs_destroy(&delete_vertex);
	igraph_destroy(&graph);

	//COMO VER EL 10% - 20% ... DE REMOSION

	return 0;
}
