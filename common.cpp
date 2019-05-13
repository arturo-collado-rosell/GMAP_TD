#ifndef COMMON
#define COMMON

#include <armadillo>
#include <complex>
#include <string>

using namespace arma;

#define PI 3.14159


void read_metaDatos(std::string file_name, int& M, double& PRF, double& sigma, double& sigmaW, int& celdas, int& I, double& lambda, int& P)
{
	FILE * pFile;
	pFile = fopen(file_name.c_str(), "r");


	float sigmaW_aux, sigma_aux, lambda_aux, PRF_aux;

	fscanf (pFile, "%i", &M);
	fscanf (pFile, "%f", &PRF_aux);
	fscanf (pFile, "%i", &celdas);
	fscanf (pFile, "%i", &I);
	fscanf (pFile, "%f", &sigmaW_aux);
	fscanf (pFile, "%f", &sigma_aux);
    fscanf (pFile, "%f", &lambda_aux);
    fscanf (pFile, "%i", &P);
    
    // Convert float into double
    PRF = double(PRF_aux);
    sigma = double(sigma_aux);
    sigmaW = double(sigmaW_aux);
    lambda = double(lambda_aux);

    fclose (pFile);

}

void construyeMatriztiempos(int M, double PRF, double* tMatriz)
{
	
	double dt = 1.0/PRF;
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			tMatriz[i*M+j] = (j-i)*dt ;
		}
	}
}

cx_vec ConstruyeVector_Cplx(double* data_input, int comienzo, int fin)
{
	int M = fin-comienzo;
	
	vec real_p(M);
	vec imag_p(M);
	for (int i = 0; i < M; ++i)
	{
		real_p(i) = data_input[comienzo+i];
	}

	for (int j = 0; j < M; ++j)
	{
		imag_p(j) = data_input[comienzo+M+j];
	}

	return cx_vec(real_p,imag_p);
}

cx_mat Construye_Matrix_muchasRealizaciones(double* data_input, int M, int I)
{
	cx_mat D(M,I);

	for (int i = 0; i < I; ++i)
	{
		D.col(i) = ConstruyeVector_Cplx(data_input, 2*i*M, 2*i*M+M);
	}

	return D;
}

void agregaDatosProcesadosAlFichero(FILE* file, cx_vec x, int M)
{
	
	for (int i = 0; i < M; ++i)
	{
		double a = real(x(i));
		fwrite(&a , sizeof(double), 1, file);
	}
	for (int i = 0; i < M; ++i)
	{
		double a = imag(x(i));
		fwrite(&a, sizeof(double), 1, file);
	}
}


void AlmacenaMatrixcompleja(cx_mat A, int M, std::string name)
{
	FILE* file = fopen(name.c_str(),"w");
	for (int i = 0; i < M; ++i)
	{
			for (int j = 0; j < M; ++j)
		{
			double a = real(A(i,j));
			fwrite(&a , sizeof(double), 1, file);
		}
			for (int j = 0; j < M; ++j)
		{
			double a = imag(A(i,j));
			fwrite(&a, sizeof(double), 1, file);
		}
	}
	fclose(file);
}

vec Multiplica_Ventana(cx_vec& x, int M, char tipo)
{
	vec w(M);
	if(tipo=='B' || tipo=='b')
	{
		for (int i = 0; i < M; ++i)
		{
			w(i) = 0.42 - 0.5*std::cos(2*PI*i/(M-1)) + 0.08*std::cos(4*PI*i/(M-1));
		}
	}
	else if(tipo=='H' || tipo == 'h')
	{
		for (int i = 0; i < M; ++i)
		{
			w(i) = 0.53836 - 0.46164*std::cos(2*PI*i/(M-1));
		}
	}

	x=x%w;

	return w;
}

vec Ventana(char tipo, int M)
{
	vec w(M);
	if(tipo=='B' || tipo=='b')
	{
		for (int i = 0; i < M; ++i)
		{
			w(i) = 0.42 - 0.5*std::cos(2*PI*i/(M-1)) + 0.08*std::cos(4*PI*i/(M-1));
		}
	}
	else if(tipo=='H' || tipo == 'h')
	{
		for (int i = 0; i < M; ++i)
		{
			w(i) = 0.53836 - 0.46164*std::cos(2*PI*i/(M-1));
		}
	}
	return w;
}


	vec DEP_estimation(cx_vec x , vec w , double PRF, int M ,bool save , int I)
{
	double U = accu(w%w);
	vec spectrum(M);
	if(I==1)
	{
		x = x%w;
		spectrum = 1.0/(PRF*U) * abs(fft(x))%abs(fft(x));

		if(save)
		{
			FILE* spectrum_file = fopen("Spectrum.bin","w");
			for (int i = 0; i < M; ++i)
			{
				double a = spectrum(i);
				fwrite(&a, sizeof(double),1,spectrum_file);
			}
		
			fclose(spectrum_file);
		}
	}
	
	
	return spectrum;
}

vec DEP_estimation_matrix(cx_mat x , vec w , double PRF, int M ,bool save , int I)
{
	double U = accu(w%w);
	vec aux_sum(M);
	aux_sum.zeros();
	vec spectrum(M); 
	if(I>1)
	{
		for (int i = 0; i < I; ++i)
		{
			x.col(i) = x.col(i)%w;
			aux_sum += abs(fft(x.col(i)))%abs(fft(x.col(i))) ;
		}
		spectrum = 1.0/(PRF*U*I)*aux_sum;

	}
	else if(I==1)
	{
			
		x = x%w;
		spectrum = 1.0/(PRF*U) * abs(fft(x))%abs(fft(x));	
		
	
	}

	if(save)
	{
		FILE* spectrum_file = fopen("Spectrum.bin","w");
			for (int i = 0; i < M; ++i)
			{
				double a = spectrum(i);
				fwrite(&a, sizeof(double),1,spectrum_file);
			}
		
			fclose(spectrum_file);
	}

	return spectrum;
}

double Noise_level(vec DEP, int M, double PRF, int I)
{
	//Shifteo alrededor de la frecuencia cero
	vec DEP_shifteada(M);
	if(M%2==0)
	{
		for(int i=0 ; i<M ; i++)
		{
			DEP_shifteada[i] =DEP[(M/2+i)%M]; 
		}
	}
	else
	{
		for(int i=0 ; i<M ; i++)
		{
			DEP_shifteada[i] =DEP[((M-1)/2+1+i)%M]; 
		}
	}

	bool flag= true;
	double dF = PRF/M;
	double NN = M;
	double porciento = 0.02;
	int i=0;
	double F;
	std::vector<double> sigma2;
	std::vector<double> sigma2N;
	std::vector<double> P;
	std::vector<double> Q;
	std::vector<double> R1;
	std::vector<double> R2;

	while(flag)
	{
		F =dF*NN;
		vec f(NN);
		for(int j=0 ; j<NN ; ++j)
		{
			f[j] = -F/2.0 + j*F/NN ;
		}
		//calculo de las magnitudes
		sigma2.push_back(accu(f%f%DEP_shifteada)/accu(DEP_shifteada) - accu(f%DEP_shifteada)/accu(DEP_shifteada) * accu(f%DEP_shifteada)/accu(DEP_shifteada) );
		sigma2N.push_back(F*F/12.0);
		P.push_back(accu(DEP_shifteada)/NN);
		Q.push_back(accu(DEP_shifteada%DEP_shifteada/NN) - P[i]*P[i]);
		R1.push_back(sigma2N[i]/sigma2[i]);
		R2.push_back(P[i]*P[i]/(Q[i]*I));
		//condicion para detenerse
		if((R1[i] < 1.05 && R2[i]>0.98) || NN<=10.0 )
			flag=false;
		else
		{
			int L = NN;
			double maxDEP = max(DEP_shifteada);
			double minDEP = min(DEP_shifteada);
			double threshold = accu(DEP_shifteada)/NN + porciento*(maxDEP-minDEP);
			int count = 0;
			std::vector<int> indices;
			for (int ii = 0; ii < NN; ++ii)
			{
				if(DEP_shifteada[ii] > threshold )
				{
					indices.push_back(ii);
					++count;
				}
			}
			//Elimino del vector de la DEP los elementos con indices dado por el vector indices
			
			for (int ii = 0; ii < count; ++ii)
			{
				DEP_shifteada.shed_row(indices[ii]-ii); //se eliminan las muestras
			}
			NN -=count;
			F -=count*dF;  
		}
		++i;	
	}

double NL = accu(DEP_shifteada)/NN;
return NL;

}

double Estima_Potencia(vec DEP, int M, double sigma, double df)
{	
	double P = std::sqrt(2.0*PI)*sigma*(DEP(M-1) + DEP(0) + DEP(1) )/(1.0 + 2.0*std::exp(-df*df/(2.0*sigma*sigma))) ;
	return P;
}

mat Rc_matrix(double P, double sigma, int M, double PRF)
{
	double* t = (double*)malloc(sizeof(double*)*M*M);
	construyeMatriztiempos( M,  PRF,  t);
	mat R(M,M);
	double Ts = 1.0/PRF;
	//rowvec R(M);
	for (int i = 0; i < M; ++i)		
	{
		for (int j = 0; j < M; ++j)
		{
			R(i,j) = P*std::exp( -2.0*PI*PI*sigma*sigma*t[i*M+j]*t[i*M+j]);
		}
		
	}
	free(t);
	return  R;


}

cx_mat Rp_matrix(double Sp, double sigma_p, double fm, int M, double PRF)
{
	double Ts = 1.0/PRF;
	double* t = (double*)malloc(sizeof(double*)*M*M);
	construyeMatriztiempos( M,  PRF,  t);
	Mat<cx_double> Rp(M,M);
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			Rp(i,j) = cx_double(Sp*std::exp(-2.0*PI*PI*sigma_p*sigma_p*t[i*M+j]*t[i*M+j])* std::cos(-2.0*PI*fm*t[i*M+j]) , Sp*std::exp(-2.0*PI*PI*sigma_p*sigma_p*t[i*M+j]*t[i*M+j])* std::sin(-2.0*PI*fm*t[i*M+j]) );
		}
	}
	free(t);
	return Rp;
	
}

template<typename T>
T sgn(T x)
{
	if(x>0)
		return T(1);
	else 
		return T(-1);
}


std::vector<double> CalculaMomentosPPP(cx_mat Ry , int M, double PRF, double NoiseP)
{
	double Pf;
	double R0=0;
	for (int i = 0; i < M; ++i)
	{
		R0 = R0 + real(Ry(i,i)); 
	}

	Pf = R0/M  - NoiseP;

	std::complex<double> R1(0,0);
	double fm;
	for (int i = 0; i < M-1; ++i)
	{
		//R1 = R1 + Ry(i,i+1);
		R1 = R1 + Ry(i+1,i);
	}
	R1 = R1/(M-1.0);
	fm = PRF/(2*PI)*arg(R1);

	double sigma;
	if(Pf/abs(R1) <= 1.0)
	{
		 sigma = 0.0;
	}
	else
	{
		   sigma = PRF/(PI*std::sqrt(2.0))* std::sqrt(std::log(Pf/abs(R1))) ;
	}
	
	std::vector<double> out(3);
	out[0] = Pf; out[1]=fm; out[2]=sigma;
	return out;
}




void saveMoments(std::vector<std::vector<double>> momentos, int Lr, std::string name_file )
{
	FILE* file = fopen(name_file.c_str(), "w");
	double aux;
	for (int i = 0; i < Lr; ++i)
	{	
		for (int j = 0; j < 3; ++j)
		{
			aux = momentos[i][j];
			fwrite(&aux,sizeof(double),1,file);
		}
	}	

	fclose(file);
}






#endif