Instrucciones de ejecución del programa:

Hay 7 archivos que deben estar todos en la misma carpeta para la ejecución del programa principal.
Se adjuntan dos "formatos". Un .rar con solo los archivos a ejecutar y un .rar con la carpeta de la libreria a utilizar dentro de la cual se encontraran dichos archivos.

* Instalar igraph. El/Los computadores donde se correrá el programa deben tener instalado la libreria igraph para C.
* Esta se puede encontrar (junto con las instrucciones de instalación) en: http://igraph.org/c/

Una vez instalada la libreria se deben agregar, a la carpeta igraph-0.7.1 (en caso de no utilizar el rar que contiene la carpeta mencionada), los 6 archivos principales.

Para compilar el programa se debe utilizar el comando gcc de la siguiente manera:

gcc -O3 -std=c99 server-graph.c -I/usr/local/include/igraph -L/usr/local/lib -ligraph -lm -pthread -o server-graph

donde server-graph.c es el nombre del archivo principal encargado de ejecutar todos los códigos. Además, -I/... y -L/... pueden variar según el lugar donde se encuentre ubicado igraph. Para esto puede ejecutar el siguiente comando:

pkg-config --libs --cflags igraph

Este le mostrará tanto el -I como el -L que deberá usar para la compilación.

Una vez compilado el programa se debe ejecutar de la siguiente manera:

./server-graph threads alpha nodes edges elem-comp rad-CI

donde:

* threads: Número de threads que se desea utilizar para la paralelización.
* alpha: Es el parámetro (exponente) de power-law. Debe ser mayor a 2.
* nodes: Número de nodos que tendrá el grafo a generar.
* edges: Número de aristas que tendrá el grafo a generar.
* elem-comp: Cantidad de elementos que se desea comparar en Kendall.
* rad-CI: Radio de la vecindad de Collective Influence.

**** NOTA *****
SE HA PROBADO CON 4 THREADS Y CON GRAFOS DE 1000 NODOS. SE PROBÓ CON UNO DE 10000 Y FUNCIONAN LAS METRICAS: COREHD Y DEGREE (las otras se demoran harto así que no corroboré el buen funcionamiento con este grafo).

LOS RESULTADOS DE KENDALL ESTARÁN EN EL ARCHIVO Kendall_Correlation.csv
ESTE TENDRÁ EL NOMBRE-ARCHIVO-GRAFO, ELEMENTOS-A-COMPARAR, COEFICIENTES-SEGUN-METRICAS

LOS GRAFOS ALEATORIOS GENERADOS ESTARAN EN redTestX.edges DONDE X CORRESPONDE AL N° DEL THREAD QUE LO GENERÓ!! 

EL METODO USADO PARA ELIMINAR LOS NODOS RESTANTES DE CI FUE "order" ES DECIR, SE ELIMINAN EN ORDEN LOS NODOS QUE VAN QUEDANDO. ESTAN LAS OPCIONES DE ELIMINACION "random" Y "degree" (eliminan de manera random los nodos que quedan y elimina según el de mayor grado, respectivamente). EL PROGRAMA PRINCIPAL SOLO CONSIDERA LA OPCION "order" 
