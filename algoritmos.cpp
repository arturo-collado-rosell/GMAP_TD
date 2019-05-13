#ifndef ALGORITMOS
#define ALGORITMOS

#include <iostream>
#include "common.cpp"
#include <chrono>



std::vector<double> GMAP_TD(cx_mat x, int I, int M, double lambda, double PRF, double sigma)
{
	vec window = Ventana('b',M);
	//spectrum estimation
	vec DEP =  DEP_estimation_matrix( x ,  window ,  PRF,  M , false, I);	
	//estima potencia
	double df = PRF/M;

	double P = Estima_Potencia(DEP,  M, sigma, df);

	//calculo de la matrix de correlacion
	cx_mat R_e(M,M);
	R_e.zeros();
	for (int i = 0; i < I; ++i)
	{
		R_e +=  x.col(i) * (x.col(i).t());
	}

	R_e = R_e/double(I);

	//nivel de ruido (Por el momento le coloco el valor simulado, estoy esperando a que se defina a nivel de grupo como estimar el piso de ruido)
	double noise = 5.0;

	//construyo la matrix A (filtro)
	mat Rc= Rc_matrix( P,  sigma, M, PRF);

	bool flag_to_invert = true;

	cx_mat A(M,M);

	if(flag_to_invert)
	{
		cx_vec eigval;
		cx_mat  eigvec;
		mat X =   Rc/(noise*PRF) + eye<mat>(M,M);
		eig_gen( eigval, eigvec,X); 

		cx_mat aux = sqrtmat(diagmat(eigval));

		cx_mat usandoDiagonalizacion = eigvec * aux * inv(eigvec) ;
	
		A = inv_sympd(usandoDiagonalizacion);
	}
	else
	{
		 A = inv(sqrtmat(Rc/(noise*PRF) + eye<mat>(M,M)));
		//cx_mat A = inv_sympd(sqrtmat(Rc/(noise*PRF) + eye<mat>(M,M))); 
	}

	//salida del filtro
	cx_mat Ry = A*R_e*A.t();

	//estima los tres momentos espectrales
	std::vector<std::vector<double>> momentos(30, std::vector<double>(3)); // 30 sera el numero maximo de iteraciones permitidas


	momentos[0]= CalculaMomentosPPP(Ry ,M, PRF, noise*PRF);
	// loop de reconstruccion del espectro del fenomeno

	double del_f = 1000; //inicio esta variable con un valor exageradamente grande
	double del_Sp = 2; // idem
	cx_mat Ry_0(Ry) ; 
	int j=0;
	bool flag=true;
	//inicializo varias matrices utiles
	cx_mat Rp_aux(M,M);
	cx_mat Rpfilter(M,M);
	cx_mat dif_Rp(M,M);

	if(momentos[0][2]!=0.0 )
	{
		while((del_f > 0.005*PRF/2.0 || del_Sp > std::pow(10,0.01) ) && j<30)
		{
			Rp_aux = Rp_matrix( momentos[j][0], momentos[j][2], momentos[j][1], M, PRF);
			Rpfilter = A*Rp_aux*A.t();
			dif_Rp = Rp_aux - Rpfilter;
			Ry = Ry_0 + dif_Rp;

			if(j >=2 )
			{
				del_f = std::abs(momentos[j][1]-momentos[j-1][1]);
				del_Sp = momentos[j][0]/momentos[j-1][0];
			}

			++j;
			momentos[j] = CalculaMomentosPPP(Ry ,M, PRF, noise*PRF);
		

			if(momentos[j][2]==0)
			{
				//std::cout<<"sigma imaginario. Este esultado no se toma\n";
				flag=false;
				break;
			}
		}
	}
	else
	{
		flag=false;
	}

	
 	if(flag && j>1)
	 	return momentos[j];
 	
    else
    	return momentos[0];
			
}




std::vector<double> GMAP_TD_Sebastian(cx_mat x, int I, int M, double PRF, double sigma, cx_vec eigval, cx_mat  eigvec)
{
	
	vec window = Ventana('b',M);
	//spectrum estimation
	vec DEP =  DEP_estimation_matrix( x ,  window ,  PRF,  M , false, I);
	//estima potencia
	double df = PRF/M;

	double P = Estima_Potencia(DEP,  M, sigma, df);
	//calculo de la matrix de correlacion
	cx_mat R_e(M,M);
	R_e.zeros();
	for (int i = 0; i < I; ++i)
	{
		R_e +=  x.col(i) * (x.col(i).t());
	}

	R_e = R_e/double(I);

	
	//nivel de ruido
	//double noise = 5.0;
	double noise = Noise_level( DEP,  M,  PRF,  I);

	//construyo la matrix A (filtro)
	cx_vec D1 = cx_vec(M);
	for (int ii = 0; ii < M; ++ii)
	{
		D1[ii] = std::sqrt(1.0/(1.0+ P/(noise*PRF)*real(eigval[ii]))) ; //Simplificacion de Sebastian
	}

	cx_mat aux = diagmat(D1);
	cx_mat A = eigvec * aux * trans(eigvec) ;
	
	//salida del filtro
	cx_mat Ry = A*R_e*A.t();
	//estima los tres momentos espectrales
	std::vector<std::vector<double>> momentos(30, std::vector<double>(3));


	momentos[0]= CalculaMomentosPPP(Ry ,M, PRF, noise*PRF);

	// loop de reconstruccion del espectro del fenomeno
	double del_f = 1000;
	double del_Sp = 2;
	cx_mat Ry_0(Ry) ; 
	int j=0;
	bool flag=true;

	cx_mat Rp_aux(M,M);
	cx_mat Rpfilter(M,M);
	cx_mat dif_Rp(M,M);

	if(momentos[0][2]!=0.0 )
	{
		while((del_f > 0.005*PRF/2.0 || del_Sp > std::pow(10,0.01) ) && j<30)
		{
			momentos[j]= CalculaMomentosPPP(Ry ,M, PRF, noise*PRF);
			//std::cout<< momentos[j][0]<<"  "<<momentos[j][1]<<"  "<<momentos[j][2]<<"\n";
			Rp_aux = Rp_matrix( momentos[j][0], momentos[j][2], momentos[j][1], M, PRF);
			Rpfilter = A*Rp_aux*A.t();
			dif_Rp = Rp_aux - Rpfilter;
			Ry = Ry_0 + dif_Rp;
			//agrego esta linea
			//Ry = Ry + dif_Rp;
			//Ry = A*Ry*A.t(); 

			if(j >=1 )
			{
				del_f = std::abs(momentos[j][1]-momentos[j-1][1]);
				del_Sp = std::abs(momentos[j][0]/momentos[j-1][0]);
			}			
			//std::cout<<"sigma ="<<momentos[j][2]<<"\n";
			if(momentos[j][2]==0)
			{
				
				flag=false;
				break;
			}
			
			//std::cout<<"del_f "<<del_f<<" del_Sp "<< del_Sp<<" \n";
			++j;	
		}
	}
	else
	{
		flag=false;
	}

 	if(flag && j>1)
 	{
 		//std::cout<<"ok"<<"numero de iteraciones ="<<j-1<<"\n";
	 	return momentos[j-1];
 	}
 		
 	
    else
    {
    	//std::cout<<"not ok"<<"numero de iteraciones ="<<j-1<<"\n";
    	return momentos[0];
    }
    	
	
	std::cout<<endl;
	
		
}





std::vector<double> GMAP_TD_cpp(cx_mat x, int I, int M, double PRF, double sigmaW, double sigma)
{

	//Calculo de los autovalores y autovectores de la matriz de autocorrelacion del clutter
	mat Rc= Rc_matrix( 1.0,  sigma, M, PRF);
	vec eigval_real;
	mat  eigvec_real;
	//calculo los autovalores y autovectores de Rc (esto se hace aca, para poder construir el filtro A de forma mas simple)
	eig_sym( eigval_real, eigvec_real,Rc);
	

	cx_vec eigval =cx_vec(eigval_real, vec(M,fill::zeros)); 
	cx_mat eigvec =cx_mat(eigvec_real, mat(M,M,fill::zeros));

	
	vec window = Ventana('b',M);
	//spectrum estimation
	vec DEP =  DEP_estimation_matrix( x ,  window ,  PRF,  M , false, I);
	//estima potencia
	double df = PRF/M;

	double P = Estima_Potencia(DEP,  M, sigmaW, df);
	//calculo de la matrix de correlacion
	cx_mat R_e(M,M);
	R_e.zeros();
	for (int i = 0; i < I; ++i)
	{
		R_e +=  x.col(i) * (x.col(i).t());
	}

	R_e = R_e/double(I);

	
	//nivel de ruido
	//double noise = 5.0;
	double noise = Noise_level( DEP,  M,  PRF,  I);

	//construyo la matrix A (filtro)
	cx_vec D1 = cx_vec(M);
	for (int ii = 0; ii < M; ++ii)
	{
		D1[ii] = std::sqrt(1.0/(1.0+ P/(noise*PRF)*real(eigval[ii]))) ; //Simplificacion de Sebastian
	}

	cx_mat aux = diagmat(D1);
	cx_mat A = eigvec * aux * trans(eigvec) ;
	
	//salida del filtro
	cx_mat Ry = A*R_e*A.t();
	//estima los tres momentos espectrales
	std::vector<std::vector<double>> momentos(30, std::vector<double>(3));


	momentos[0]= CalculaMomentosPPP(Ry ,M, PRF, noise*PRF);

	// loop de reconstruccion del espectro del fenomeno
	double del_f = 1000;
	double del_Sp = 2;
	cx_mat Ry_0(Ry) ; 
	int j=0;
	bool flag=true;

	cx_mat Rp_aux(M,M);
	cx_mat Rpfilter(M,M);
	cx_mat dif_Rp(M,M);

	if(momentos[0][2]!=0.0 )
	{
		while((del_f > 0.005*PRF/2.0 || del_Sp > std::pow(10,0.01) ) && j<30)
		{
			momentos[j]= CalculaMomentosPPP(Ry ,M, PRF, noise*PRF);
			//std::cout<< momentos[j][0]<<"  "<<momentos[j][1]<<"  "<<momentos[j][2]<<"\n";
			Rp_aux = Rp_matrix( momentos[j][0], momentos[j][2], momentos[j][1], M, PRF);
			Rpfilter = A*Rp_aux*A.t();
			dif_Rp = Rp_aux - Rpfilter;
			Ry = Ry_0 + dif_Rp;
			//agrego esta linea
			//Ry = Ry + dif_Rp;
			//Ry = A*Ry*A.t(); 

			if(j >=1 )
			{
				del_f = std::abs(momentos[j][1]-momentos[j-1][1]);
				del_Sp = std::abs(momentos[j][0]/momentos[j-1][0]);
			}			
			//std::cout<<"sigma ="<<momentos[j][2]<<"\n";
			if(momentos[j][2]==0)
			{
				
				flag=false;
				break;
			}
			
			//std::cout<<"del_f "<<del_f<<" del_Sp "<< del_Sp<<" \n";
			++j;	
		}
	}
	else
	{
		flag=false;
	}

 	if(flag && j>1)
 	{
 		//std::cout<<"ok"<<"numero de iteraciones ="<<j-1<<"\n";
	 	return momentos[j-1];
 	}
 		
 	
    else
    {
    	//std::cout<<"not ok"<<"numero de iteraciones ="<<j-1<<"\n";
    	return momentos[0];
    }
    	
	
	std::cout<<endl;
	
		
}


#endif
