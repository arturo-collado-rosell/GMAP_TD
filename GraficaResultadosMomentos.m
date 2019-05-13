%Script para levantar los resultados obtenidos con GMAP-TD.
% close all;
clear;
I = 1000;
P = 20;
Ts = 500e-6;                   % tiempo de muestreo (PRP)
fs = 1/Ts;                      % frecuencia de muestreo
fc = 5.0e9;                     % frecuencia de portadora [GHz]
c= 3e8;
lambda = c/fc;                  % longitud de onda

vs = lambda*fs/2;               % velocidad de muestreo
vm = linspace(0,0.8*vs/2,P);

Sp = 10^(50/10);
sigma_pv = 3; 
for q=1:P
    
    file_name = strcat('archivos_out/momentos' , num2str(q-1),'.bin');
    file = fopen(file_name);
    momentos = zeros(3,I);
    potencia=0;
    velocidadm=0;
    ancho=0;
    j=1;
    
    
    for i=1:I
        momentos(:,i) = fread(file,3,'double');
        if momentos(3,i)==0
            
        else
            potencia(j)= 10*log10(momentos(1,i));
            velocidadm(j)= -momentos(2,i)*lambda/2;
            ancho(j)= momentos(3,i)*lambda/2;
            j = j+1;
            
        end
    end
    fclose(file);
    
    
    PotenciaM(q) = mean(potencia) ; stdPotenciaM(q) = std(potencia);
    VelocidadM(q) = mean(velocidadm); stdVelocidadM(q) = std(velocidadm);
    anchoE(q) = mean(ancho);  stdanchoE(q) = std(ancho);
    
    
    
    
    
    
    
%     figure;
%     subplot(3,1,1); plot(potencia);
%     subplot(3,1,2); plot(velocidadm);
%     subplot(3,1,3); plot(ancho);
    
    
    
    
end

figure; 
subplot(3,1,1) ; errorbar(vm/(vs/2), PotenciaM-10*log10(Sp), stdPotenciaM , 'LineWidth',2); ylabel('[dB]'); title('Sesgo potencia') 
subplot(3,1,2) ; errorbar(vm/(vs/2), VelocidadM-vm, stdVelocidadM, 'LineWidth',2); ylabel('m/s'); title('Sesgo velocidad media');
subplot(3,1,3) ; errorbar(vm/(vs/2), anchoE-sigma_pv, stdanchoE, 'LineWidth',2); ylabel('m/s'); xlabel('v/v_{max}') ; title('Sesgo ancho espectral');

%%
% %levanta datos de la matrix to print
% M=64;
% file1 = fopen('A1.bin');
% file2 = fopen('A2.bin');
% Matrix = zeros(M);
% for i=1:M
%     Re1 = fread(file1,M,'double');
%     Im1 = fread(file1,M,'double');
%     Matrix1(:,i) = Re1+1j*Im1;
%     
%     Re2 = fread(file2,M,'double');
%     Im2 = fread(file2,M,'double');
%     Matrix2(:,i) = Re2+1j*Im2;
%     
% end
% figure; imagesc(abs(Matrix2-Matrix1));
% 
% %%
% %estimacio ndel nivel de ruido usando la rutina de matlab
