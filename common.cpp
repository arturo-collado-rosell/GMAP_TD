#ifndef COMMON
#define COMMON

#include <armadillo>
#include <complex>
#include <string>

using namespace arma;

#define PI 3.14159

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


double Estima_Potencia(vec DEP, int M, double sigma, double df)
{	
	double P = std::sqrt(2.0*PI)*sigma*(DEP(M-1) + DEP(0) + DEP(1) )/(1.0 + 2.0*std::exp(-df*df/(2.0*sigma*sigma))) ;
	return P;
}



mat Rc_matrix(double P, double sigma, int M, double PRF)
{
	double Ts = 1.0/PRF;
	rowvec R(M);
	for (int i = 0; i < M; ++i)
	{
		R(i) = P*std::exp( -2.0*PI*PI*sigma*sigma*i*i*Ts*Ts);
	}
	return  toeplitz(R);


}

cx_mat Rp_matrix(double Sp, double sigma_p, double fm, int M, double PRF)
{
	double Ts = 1.0/PRF;
	Row<cx_double> Rp(M);
	for (int j = 0; j < M; ++j)
	{
		Rp(j) = cx_double(Sp*std::exp(-2.0*PI*PI*sigma_p*sigma_p*j*j*Ts*Ts)* std::cos(-2.0*PI*fm*j*Ts) , Sp*std::exp(-2.0*PI*PI*sigma_p*sigma_p*j*j*Ts*Ts)* std::sin(-2.0*PI*fm*j*Ts) );
		
	}

	Row<cx_double> Rpaux(M);
	Rpaux(0) = Rp(0); Rpaux(span(1,M-1)) = conj(Rp(span(1,M-1)));	
	

	return toeplitz(Rpaux,Rp);
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
	fm = -PRF/(2*PI)*arg(R1);

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




void saveMoments(std::vector<std::vector<double>> momentos, int Lr )
{
	FILE* file = fopen("momentos.bin", "w");
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