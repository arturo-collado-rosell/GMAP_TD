function [NL, umbral ] = Noiselevel(DEP,PRF,I)
%Parámetros de entrada:

% DEP es la densidad espectral de potencia
% PRF, frecuencia de muestreo
% I, número de realizaciones

%Parámetros de salida:
% NL, piso de ruido estimado
% umbral, último valor del threshold antes de calcular NL
DEP = fftshift(DEP);
flag =1 ; %condición para que se ejecute el algoritmo
N = length(DEP);
dF = PRF/N;
i=1 ;
NN = N ;
porciento = 0.05;

while flag==1
    
    F = dF*NN;
    f = -F/2:F/NN:F/2-F/NN ;
    %cálculo de los parámetros
    sigma2(i) = sum(f.^2.*DEP)/sum(DEP) - (sum(f.*DEP)/(sum(DEP)))^2 ;
    sigma2N(i) = F^2/12 ;
    P(i) = sum(DEP)/NN ;
    Q(i) = sum(DEP.^2/NN) - P(i)^2 ;
    R1(i) = sigma2N(i)/sigma2(i);
    R2(i) = P(i)^2/(Q(i)*I) ;
    
    %condición de detener el algoritmo
    if  R1(i)<1.1 && R2(i)>0.97  || (flag==0 || NN<5 )
        flag = 0 ;
    else
        L = length(DEP);
        maxDEP = max(DEP);
        minDEP = min(DEP);
        threshold = sum(DEP)/NN + porciento*(maxDEP-minDEP) ;
        p=1;
        indices = [];
        if NN==5
           flag=0; 
        end
        for j=1:L
            if DEP(j)>threshold           
                indices(p) = j;
                p = p +1;
                NN = NN-1 ;
                F = F- dF ;             
            end
        end
        DEP(indices) = [];        
    end
    i = i+1;
end

NL = sum(DEP)/NN;
umbral = threshold;


