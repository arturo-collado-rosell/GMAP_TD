//GMAP-TD implementation in c++ using armadillo library for matrix operations
//author: Arturo Collado Rosell 11/1/19

#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <armadillo>
#include <complex>
#include <chrono>

#include "common.cpp" //fichero cpp donde se encuentran funciones utiles que no forman parte de GMAP-TD
#include "algoritmos.cpp" //Aca van los algoritmos , por el momento GMAP-TD



int main(int argc, char const *argv[])
{
	
/*************** En esta versión los parámetros que se conocen se cargan manualmente. en próximas versiones esto será externo al código******/
	int M=64; //Número de muetras 
	int I=1;  // Número de realizaciones para estimar el espectro. (en la práctica I=1), con datos simulados puede ser > 1
	int celdas = 1000;  //Número de celdas en rango 
	double PRF = 2000.0; //[Hz] Frecuencia de repetición de pulsos
	double sigma_c = 23.75; // [Hz] ancho de la DEP del clutter. hay que tener en cuenta el efecto de la ventana. Este valor cambia según el número de muestras 
	double lambda = 0.06; // longitud de onda [m]

	auto start = std::chrono::system_clock::now(); // para medir tiempo
	double* buffer;
	//allocate memory to contain the whole file
	buffer = (double*)malloc(sizeof(double)*2*M*I*celdas); //alloco memoria para almacenar los datos radar
	if(buffer==NULL){fputs("Memory error", stderr);}

	//Open data file
	FILE* pFile;
	std::string name_file_to_open = "archivos_in/data" + std::to_string(0) + ".bin";
	pFile = fopen(name_file_to_open.c_str(), "r");
	pFile = fopen("data1.bin", "r");
	if(pFile==NULL) {fputs("File error", stderr);}

	//copy the file into the buffer
	size_t result = fread(buffer,sizeof(double),2*M*I*celdas, pFile);
	// Función que construye la matriz de datos IQ de Mxceldas
	cx_mat x =Construye_Matrix_muchasRealizaciones(buffer,  M,  I*celdas); // cx_mat es una forma de definir matrices usando armadillo
	
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> diff = end-start;
	std::cout << "Tiempo trascurrido para levantar los datos  y armar la matriz IQ " << diff.count() << " s\n"; // Se imprime el tiempo de levantar los datos y contruir la matriz IQ
	/*****Se termina de levantar los datos y se obtiene una matriz tal y como se estuviese trabajando en Matlab*****/

	/*******************************GMAP-TD************************************/
	//Comienzo del testeo del algoritmo GMAP-TD
	 start = std::chrono::system_clock::now(); //Inicio del tiempo
	
	//Se crea la matriz de dimensiones 3xceldas, donde se almacenarán los resultados de los 3 primeros momentos del fenómeno
	std::vector<std::vector<double>> resultado(celdas, std::vector<double>(3));
	//Para cada celda se realiza el llamado a la rutina GMAP_TD
	for (int i = 0; i < celdas; ++i)
	{
		resultado[i]=GMAP_TD(x.col(i),I,M, lambda, PRF, sigma_c);
		//GMAP_TD(x.col(i),I,M, lambda, PRF, sigma_c);
		//std::cout<<i<<std::endl;
	}

	 end = std::chrono::system_clock::now();
	 diff = end-start;
	std::cout << "Tiempo trascurrido usando chrono " << diff.count() << " s\n";
	
	//save moments into a binary file 
	std::string name_file = "archivos_out/momentos" + std::to_string(0) + ".bin"; 
	saveMoments(resultado,celdas,name_file);

	free(buffer);
	return 0;
}