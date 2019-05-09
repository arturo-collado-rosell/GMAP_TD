



#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <armadillo>
#include <complex>
#include <chrono>

#include "common.cpp" //fichero cpp donde se encuentran funciones utiles que no forman parte de GMAP-TD
#include "algoritmos.cpp" //Aca van los algoritmos , por el momento GMAP-TD


int M=64; //Número de muetras 
int I=1;  // Número de realizaciones para estimar el espectro. (en la práctica I=1), con datos simulados puede ser > 1
int celdas = 1000;  //Número de celdas en rango 
double PRF = 2000.0; //[Hz] Frecuencia de repetición de pulsos
double sigma_c = 23; // [Hz] ancho de la DEP del clutter. hay que tener en cuenta el efecto de la ventana. Este valor cambia según el número de muestras 
double sigma = 8.33;    //[Hz] el ancho teorico del clutter
double lambda = 0.06; // longitud de onda [m]
int P = 20; //numero de velocidades distintas



int main(int argc, char const *argv[])
{
	
	for (int i = 0; i < P; ++i)
	{
		
		double* buffer;
		//allocate memory to contain the whole file
		buffer = (double*)malloc(sizeof(double)*2*M*I*celdas); //alloco memoria para almacenar los datos radar
		if(buffer==NULL){fputs("Memory error", stderr);}

		//Open data file
		FILE* pFile;
		std::string name_file_to_open = "archivos_in/data" + std::to_string(i+1) + ".bin";
		pFile = fopen(name_file_to_open.c_str(), "r");
		if(pFile==NULL) {fputs("File error", stderr);}

		//copy the file into the buffer
		size_t result = fread(buffer,sizeof(double),2*M*I*celdas, pFile);
		// Función que construye la matriz de datos IQ de Mxceldas
		cx_mat x =Construye_Matrix_muchasRealizaciones(buffer,  M,  I*celdas); // cx_mat es una forma de definir matrices usando armadillo
	
	
		/*****Se termina de levantar los datos y se obtiene una matriz tal y como se estuviese trabajando en Matlab*****/

		/*******************************GMAP-TD************************************/
		//Comienzo del testeo del algoritmo GMAP-TD
		auto start = std::chrono::system_clock::now(); //Inicio del tiempo
	
		//Se crea la matriz de dimensiones 3xceldas, donde se almacenarán los resultados de los 3 primeros momentos del fenómeno
		std::vector<std::vector<double>> resultado(celdas, std::vector<double>(3));
	
		//Construcción de la matriz Rc Sin tener en cuenta la potencia del clutter.
		mat Rc= Rc_matrix( 1.0,  sigma, M, PRF);
		vec eigval_real;
		mat  eigvec_real;
		//calculo los autovalores y autovectores de Rc (esto se hace aca, para poder construir el filtro A de forma mas simple)
		eig_sym( eigval_real, eigvec_real,Rc);
	

		cx_vec eigval =cx_vec(eigval_real, vec(M,fill::zeros)); 
		cx_mat eigvec =cx_mat(eigvec_real, mat(M,M,fill::zeros)); 

		//Para cada celda se realiza el llamado a la rutina GMAP_TD
		for (int i = 0; i < celdas; ++i)
		{
			//std::cout<<i<<"\n";
			if(I==1)
				resultado[i]=GMAP_TD_Sebastian(x.col(i),I,M, lambda, PRF, sigma_c, eigval, eigvec);
			else
				resultado[i]=GMAP_TD_Sebastian(x,I,M, lambda, PRF, sigma_c, eigval, eigvec);
		}

		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> diff = end-start;
		std::cout << "Tiempo trascurrido usando chrono " << diff.count() << " s\n";
	
		//save moments into a binary file
		std::string name_file = "archivos_out/momentos" + std::to_string(i) + ".bin"; 
		saveMoments(resultado,celdas,name_file);

		free(buffer);
		
	}
return 0;
}