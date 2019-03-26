%Script para levantar los resultados obtenidos con GMAP-TD.
% close all;
lambda=0.06;
I = 100;
file = fopen('momentos.bin');
momentos = zeros(3,I);
potencia=0;
velocidadm=0;
ancho=0;
j=1;


for i=1:I
    momentos(:,i) = fread(file,3,'double');
    if momentos(3,i)==0
        
    else
        potencia(j)= momentos(1,i);
        velocidadm(j)= momentos(2,i);
        ancho(j)= momentos(3,i);
        j = j+1;
        
    end
end

fclose(file);

figure; plot(potencia);
figure; plot(velocidadm*lambda/2);
figure; plot(ancho*lambda/2);


%%
%levanta datos de la matrix to print
M=64;
file1 = fopen('A1.bin');
file2 = fopen('A2.bin');
Matrix = zeros(M);
for i=1:M
    Re1 = fread(file1,M,'double');
    Im1 = fread(file1,M,'double');
    Matrix1(:,i) = Re1+1j*Im1;
    
    Re2 = fread(file2,M,'double');
    Im2 = fread(file2,M,'double');
    Matrix2(:,i) = Re2+1j*Im2;
    
end
figure; imagesc(abs(Matrix2-Matrix1));