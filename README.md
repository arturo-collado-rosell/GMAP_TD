1- Este repo fue creado con la intención de implementar GMAP-TD en c++ (paso previo a su implementación en cuda), con el objetivo de encontrar cuellos de botellas en la ejecución etc.

2- Una particularidad de la implementación es que se usó la librería armadillo, la cual es una librería especializada en algebra lineal y es lo mas similar a Matlab que se encuentra en c++. Lo único que por el momento es imprescindible de esta libreria es la función para calcular autovalores y autovectores de una matriz. El resto de las cosas se pueden implementar a mano (multiplicacion de matrices, raiz cuadrada matricial, etc).

3- La estructura de los archivos es la siguiente. Por un lado tenemos los archivos .m, los cuales son archivos Matlab cuyo objetivo es simular los datos sintéticos y también levantar los resultados provenientes del procesamiento y graficarlos. Por otro lado, están los archivos .cpp, en los cuales está la implementación y testeo de GMAP-TD. A continuación el nombre del archivo con su breve descripción.

   1- genera_datos.m -> genera los datos sintéticos. Utiliza el repo LibraryMeteo (addpath('../LibraryMeteo/dataGen'); addpath('../LibraryMeteo/spectrumEstimate');)

   2- GMAP_TD_test.cpp -> función para testear la implementación de GMAP-TD usando la construcción del filtro de forma naive (sin la nueva implementación de Sebastian).

   3- GMAP_TD_test_modificacionSebastian.cpp ->  función para testear la implementación de GMAP-TD usando la construcción del filtro de forma eficiente, con la modificación que hizo Sebastian.

   4- algoritmos.cpp -> aca se encuentra la implementación del algoritmo GMAP-TD. Se encuentras las dos posibles implementaciones, usando la modificación de sebastian y la que usa la forma naive(calculando la inversa de la raiz matricial ... , etc). Si bien, la nueva implementación es mucho más rápida y es la que se usará, deje las dos.  

   5-common.cpp -> es una archivo que contiene varias funciones útiles y necesarias, las cuales son invocadas por el algoritmo GMAP-TD.

   6- GraficaResultadosMomentos.m -> script de matlab que levanta el resultado de los momentos estimados por GMAP-TD y los grafica.

Forma de Almacenar los datos simulados.
 Los datos simulados en matlab se almacenan en forma binaria para que sea más facil levantar los datos en c++.	
 
 Ya esta implementada una version de la estimacion de ruido usando el paper del 1973 Objective Determination of the Noise Level in Doppler Spectra
 
 Como compilar
    Para compílar necesitan tener instalado la librería armadillo, si uno busca en http://arma.sourceforge.net/ encontrará la información necesaria para instalarla. Yo uso linux y compilo de la siguiente manera
    en la consola.
    g++ -std=c++11 nombre_archivo_a_compilar.cpp -o nombre_archivo_de_salida -larmadillo
