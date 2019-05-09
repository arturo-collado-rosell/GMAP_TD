%% Script de generacion de datos simulados para ser procesados por GMAP_TD inplementado en c++
%% 25/03/2019. Arturo
%% Generación de datos

clc;clear;close all;

addpath('../LibraryMeteo/dataGen');
addpath('../LibraryMeteo/spectrumEstimate');

situacion = 'fenomenoclutteryruido';

situacionposible = {'fenomeno',...
    'fenomenoyruido',...
    'fenomenoclutteryruido'};

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Definicion de parametros
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

Ts = 500e-6;                   % tiempo de muestreo (PRP)
fs = 1/Ts;                      % frecuencia de muestreo

fc = 5.0e9;                     % frecuencia de portadora [GHz]
c = 3e8;                        % velocidad de la luz en el vacio [m/s];

lambda = c/fc;                  % longitud de onda

vs = lambda*fs/2;               % velocidad de muestreo

CSR = 40;                       % relación clutter a señal
SNR = 20;                       % relación señal a ruido

Sp = 10^(50/10);                         % potencia del fenomeno (mantener en 1)
% vm = 0.05*vs;                   % velocidad media del fenomeno [m/s]
vm = linspace(0,0.8*vs/2,20);
sigma_pv = 3;                   % ancho espectral del fenomeno [m/s]


sigma_pf = 2*sigma_pv/lambda;   % ancho espectral del fenomeno [Hz]
var_pf = sigma_pf.^2;           % varianza del espectro del fenomeno

Sc = (10^(0.1*CSR))*Sp;         % potencia del clutter
sigma_cv = 0.25;                % ancho espectral del clutter [m/s]

sigma_cf = 2*sigma_cv/lambda;   % ancho espectral del clutter [Hz]
var_cf = sigma_cf.^2;           % varianza del espectro del clutter

No2 = Sp*Ts/(10^(0.1*SNR));     % No/2

M = 64;                         % nro de muestras de la señal
I = 1000;                          % nro de realizaciones a generar

ventana = 'blackman';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Carga de parametros
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
for qq =1:length(vm)
    fm = -2*vm(qq)/lambda;
    switch situacion
        
        case situacionposible{1}
            
            signalParams.phenom.pwr = Sp;
            signalParams.phenom.vme = vm(qq);
            signalParams.phenom.vva = sigma_pv.^2;
            
            signalParams.lambda = lambda;
            
            genParams.nSamples = M;
            genParams.nRealiz = I;
            genParams.fs = fs;
            %         genDefault.samplesFactor
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
            % Valores analiticos para comparacion
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
            
            % forma del espectro deseado para el grafico
            % frecuencia media del fenomeno [Hz]
            MM = 10*M;
            ff = (-(MM/2):(MM/2-1))*fs/MM;
            S_f = theorSpecSampled(Sp,fm,var_pf,ff,[-2:2]);
            
            S_v = 2/lambda*S_f;
            vv = -lambda*ff/2;
            
            Fdep = figure;P = plot(vv/vs,10*log10(S_v),'k');hold on;
            set(P,'LineWidth',2);
            
            % autocorrelacion deseada
            m = 0:(M-1);
            R = Sp*exp(-2*(pi^2)*var_pf*((m*Ts).^2)).*exp(-1i*2*pi*fm*m*Ts);
            
            tau = ((-M+1):(M-1))*Ts;
            
            Fcorr = figure;P = plot(tau/Ts,abs([R(end:-1:2) R]),'r');hold on;
            set(P,'LineWidth',2);
            
            % calculo de la media de la estimacion del espectro
            spectParams.sampleTime = Ts;
            spectParams.window = ventana;
            spectParams.correlation = R;
            
            % Como me interesa solo la media, paso en datos un vector de ceros.
            specDen = powerSpectDenEst(zeros(1,M),spectParams);
            
            media_peridiograma_v = 2/lambda * specDen.depEstMean;
            v = -lambda * specDen.freq /2;
            
            figure(Fdep); P = plot(v/vs,10*log10(media_peridiograma_v),'r');
            set(P,'LineWidth',2);
            
            % pdf de primer orden deseada
            variq = Sp/2;
            
            dz = 0.01;
            zb = (-3:dz:3)*sqrt(variq);
            
            fz = 1/sqrt(2*pi*variq)*exp(-(zb.^2)/(2*variq));
            
            Fpdf = figure;
            subplot(1,2,1);P = plot(zb,fz,'r');hold on;set(P,'LineWidth',2);
            subplot(1,2,2);P = plot(zb,fz,'r');hold on;set(P,'LineWidth',2);
            
        case situacionposible{2}
            
            signalParams.phenom.pwr = Sp;
            signalParams.phenom.vme = vm(qq);
            signalParams.phenom.vva = sigma_pv.^2;
            
            signalParams.noise.dep = No2;
            
            signalParams.lambda = lambda;
            
            genParams.nSamples = M;
            genParams.nRealiz = I;
            genParams.fs = fs;
            %         genDefault.samplesFactor
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
            % Valores analiticos para comparacion
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
            
            % forma del espectro deseado para el grafico
            MM = 10*M;
            ff = (-(MM/2):(MM/2-1))*fs/MM;
            S_f = theorSpecSampled(Sp,fm,var_pf,ff,[-2:2])  + No2;
            
            S_v = 2/lambda*S_f;
            vv = -lambda*ff/2;
            
            Fdep = figure;P = plot(vv/vs,10*log10(S_v),'k');hold on;
            set(P,'LineWidth',2);
            
            % autocorrelacion deseada
            m = 0:(M-1);
            R = Sp*exp(-2*(pi^2)*var_pf*((m*Ts).^2)).*exp(-1i*2*pi*fm*m*Ts) + (No2/Ts)*(m==0);
            
            tau = ((-M+1):(M-1))*Ts;
            
            Fcorr = figure;P = plot(tau/Ts,abs([R(end:-1:2) R]),'r');hold on;
            set(P,'LineWidth',2);
            
            % calculo de la media de la estimacion del espectro
            spectParams.sampleTime = Ts;
            spectParams.window = ventana;
            spectParams.correlation = R;
            
            % Como me interesa solo la media, paso en datos un vector de ceros.
            specDen = powerSpectDenEst(zeros(1,M),spectParams);
            
            media_peridiograma_v = 2/lambda * specDen.depEstMean;
            v = -lambda * specDen.freq /2;
            
            figure(Fdep); P = plot(v/vs,10*log10(media_peridiograma_v),'r');
            set(P,'LineWidth',2);
            
            % pdf de primer orden deseada
            variq = (Sp + No2/Ts)/2;
            
            dz = 0.01;
            zb = (-3:dz:3)*sqrt(variq);
            
            fz = 1/sqrt(2*pi*variq)*exp(-(zb.^2)/(2*variq));
            
            Fpdf = figure;
            subplot(1,2,1);P = plot(zb,fz,'r');hold on;set(P,'LineWidth',2);
            subplot(1,2,2);P = plot(zb,fz,'r');hold on;set(P,'LineWidth',2);
            
            
        case situacionposible{3}
            
            signalParams.phenom.pwr = Sp;
            signalParams.phenom.vme = vm(qq);
            signalParams.phenom.vva = sigma_pv.^2;
            
            signalParams.noise.dep = No2;
            
            signalParams.clutter.pwr = Sc;
            signalParams.clutter.vva = sigma_cv.^2;
            
            signalParams.lambda = lambda;
            
            genParams.nSamples = M;
            genParams.nRealiz = I;
            genParams.fs = fs;
            %         genDefault.samplesFactor
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
            % Valores analiticos para comparacion
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
            
            % forma del espectro deseado para el grafico
            MM = 10*M;
            ff = (-(MM/2):(MM/2-1))*fs/MM;
            
            S_f = theorSpecSampled(Sc,0,var_cf,ff,[-2:2])  + ...
                theorSpecSampled(Sp,fm,var_pf,ff,[-2:2]) + ...
                No2;
            
            S_v = 2/lambda*S_f;
            vv = -lambda*ff/2;
            
            Fdep = figure;P = plot(vv/vs,10*log10(S_v),'k');hold on;
            set(P,'LineWidth',2);
            
            % autocorrelacion deseada
            m = 0:(M-1);
            R = Sc*exp(-2*(pi^2)*var_cf*((m*Ts).^2)) + ...
                Sp*exp(-2*(pi^2)*var_pf*((m*Ts).^2)).*exp(-1i*2*pi*fm*m*Ts) + (No2/Ts)*(m==0);
            
            tau = ((-M+1):(M-1))*Ts;
            
            Fcorr = figure;P = plot(tau/Ts,abs([R(end:-1:2) R]),'r');hold on;
            set(P,'LineWidth',2);
            
            % calculo de la media de la estimacion del espectro
            spectParams.sampleTime = Ts;
            spectParams.window = ventana;
            spectParams.correlation = R;
            
            % Como me interesa solo la media, paso en datos un vector de ceros.
            specDen = powerSpectDenEst(zeros(1,M),spectParams);
            
            media_peridiograma_v = 2/lambda * specDen.depEstMean;
            v = -lambda * specDen.freq /2;
            
            figure(Fdep); P = plot(v/vs,10*log10(media_peridiograma_v),'r');
            set(P,'LineWidth',2);
            
            % pdf de primer orden deseada
            variq = (Sc + Sp + No2/Ts)/2;
            
            dz = 0.01;
            zb = (-3:dz:3)*sqrt(variq);
            
            fz = 1/sqrt(2*pi*variq)*exp(-(zb.^2)/(2*variq));
            
            Fpdf = figure;
            subplot(1,2,1);P = plot(zb,fz,'r');hold on;set(P,'LineWidth',2);
            subplot(1,2,2);P = plot(zb,fz,'r');hold on;set(P,'LineWidth',2);
            
        otherwise
            
            error('No existe está situación');
            
    end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    % Generacion de los datos
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    
    z = dataGenWeatherlikeDopplerSpectra(signalParams,genParams);
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    % Estimación de la DEP de los datos generados
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    
    spectParams.sampleTime = Ts;
    spectParams.window = ventana;
    specDen = powerSpectDenEst(z,spectParams);
    
    Szz_v = 2/lambda*specDen.depEst;
    v = -lambda*specDen.freq/2;
    
    figure(Fdep);P = plot(v/vs,10*log10(Szz_v),'b');grid;
    set(P,'LineWidth',2);
    eje = axis; axis([eje(1) eje(2) -30 eje(4)]);
    xlabel('Velocidad normalizada (v/vs)');ylabel('[dB]');
    legend('DEP analitica','Media estimador DEP','DEP estimada','Location','Best');
    
    
    %%Noise level
    [NL, umbral ] = Noiselevel(specDen.depEst,fs,I)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Save data into a binary file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if I==1
        
        fileID = fopen(strcat('data', num2str(qq), '.bin'),'w');
        
        fwrite(fileID, real(z),'double');
        fwrite(fileID, imag(z),'double');
        
        fclose(fileID);
        
    else
        fileID = fopen(strcat('archivos_in/data', num2str(qq), '.bin'),'w');
        for i=1:I
            fwrite(fileID, real(z(i,:)),'double');
            fwrite(fileID, imag(z(i,:)),'double');
        end
        
        fclose(fileID);
        
    end
    
    
end

