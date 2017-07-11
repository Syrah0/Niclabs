#include "igraph.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void free_complist(igraph_vector_ptr_t *complist) {
  long int i;
  for (i=0; i<igraph_vector_ptr_size(complist); i++) {
    igraph_destroy(VECTOR(*complist)[i]);
    free(VECTOR(*complist)[i]);
  }
}

int main()
{
int i,j,k,l;

FILE *F;
char filename[32];

igraph_t g,gaux;
igraph_vector_t eb,google,google2, fb, values, csize; 
igraph_matrix_t laplace;

igraph_vector_init(&google, 1);
igraph_vector_init(&google2, 1);
igraph_vector_init(&fb, 1);
igraph_vector_init(&csize, 1);
igraph_vector_init(&values, 800);

igraph_matrix_init(&laplace,1,1);

F = fopen("red3.edges","r");
igraph_read_graph_edgelist(&g, F, 0, 0);
fclose (F);

igraph_simplify(&g, 1, 1, 0);
igraph_vector_init(&eb, igraph_ecount(&g));

fprintf(stderr,"Start...\n");

j = 0;
l= 1000;
while (igraph_vector_size(&eb) > 0)
  {
  double avg=0.0;
  
/*  Subcomponente mas grande */
  igraph_clusters(&g, NULL, &csize, NULL, IGRAPH_WEAK);
  //fprintf (stdout,"%d %1.2f ",j++,igraph_vector_max(&csize));
  
  if (igraph_vector_max(&csize) < l)
    {
    l = igraph_vector_max(&csize);
    sprintf (filename,"grafo%d.edges",j++);
    if (j % 79 == 0) 
      {
      F = fopen(filename,"w");
      igraph_write_graph_edgelist(&g,F);
      fclose(F);
      fprintf (stderr,"file: %s with %d maxCC-nodes\n",filename,l);
      }
    }
  /* Calculo de eccentricity 
  igraph_eccentricity(&g, &google, igraph_vss_1(183), IGRAPH_ALL);
  igraph_eccentricity(&g, &google2, igraph_vss_1(447), IGRAPH_ALL);
  igraph_eccentricity(&g, &fb, igraph_vss_1(411), IGRAPH_ALL);
  fprintf (stdout,"%1.2f %1.2f %1.2f ",VECTOR(google)[0],VECTOR(google2)[0],VECTOR(fb)[0]);

  /* Calculo de net connectivity 

  igraph_subcomponent(&g, &google, 183, IGRAPH_ALL);
  igraph_subcomponent(&g, &google2, 447,IGRAPH_ALL);
  igraph_subcomponent(&g, &fb, 411,IGRAPH_ALL);

  fprintf (stdout,"%ld %ld %ld ",igraph_vector_size(&google),igraph_vector_size(&google2),igraph_vector_size(&fb));
  /* Algebraic conectivity 

  igraph_induced_subgraph(&g,&gaux,igraph_vss_vector(&google),IGRAPH_SUBGRAPH_AUTO);
  igraph_laplacian(&gaux, &laplace, NULL, 0, NULL);
  igraph_lapack_dsyevr(&laplace, IGRAPH_LAPACK_DSYEV_ALL, 0, 0, 0, 0, 0, 1e-10, &values, 0, 0);
  igraph_vector_sort(&values);

  igraph_degree(&gaux,&google2,igraph_vss_all(),IGRAPH_ALL, 0);
  fprintf (stdout,"%1.4f %1.4f %1.4f\n", igraph_vector_min(&google2),((1.0*igraph_vcount(&gaux))/(igraph_vcount(&gaux)-1.0))*igraph_vector_min(&google2), VECTOR(values)[1]);
/*
  fprintf (stdout,"%1.4f %1.4f \n", VECTOR(values)[1], 1.1);
*/

  /* Buscamos el maximo */
  igraph_edge_betweenness(&g, &eb, IGRAPH_UNDIRECTED, /*weights=*/ 0);
  igraph_delete_edges(&g,igraph_ess_1(igraph_vector_which_max(&eb)));
  igraph_simplify(&g, 1, 1, 0);
  
  igraph_destroy(&gaux);
  }

igraph_vector_destroy(&eb);
igraph_destroy(&g);
return 0;
} 
