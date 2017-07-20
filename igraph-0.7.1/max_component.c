/* Se encarga de calcular la componente conexa mas grande */
igraph_t* max_component(igraph_t *g){
	igraph_vector_ptr_t complist;

	/* Se inicializa vector que contendra todas las componentes conexas del grafo */
	igraph_vector_ptr_init(&complist, 1);
	igraph_decompose(g,&complist,IGRAPH_WEAK,-1,2); // se realiza la descomposicion del grafo en sus componentes

	/* Se verifica cual es la componente mas grande */
	int max = 0;
	for(int i=1; i<igraph_vector_ptr_size(&complist); i++){
		if(igraph_vector_size(VECTOR(complist)[max]) < igraph_vector_size(VECTOR(complist)[i])){ // Hay una mas grande que la componente actual
			max = i;
		}
	}
	return VECTOR(complist)[max]; // componente conexa mas grande
}