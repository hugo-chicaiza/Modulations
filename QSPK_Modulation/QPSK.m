%% PROCESAMIENTO SEÑAL
f=1; %Frecuencia de reloj binario 
% Normalización de la señal
load('NRC12225.mat')
sam=fs; % Muestras por simbolo
t=linspace(0,2*pi,sam);


senal_normalizada = normalize(Senal); % Normalizar la señal
nivel_bits = 2; % Número de bits por símbolo en el esquema QPSK
num_simbolos = floor(length(senal_normalizada) / nivel_bits); % Calcular el número de símbolos
bits_cuantificados = zeros(num_simbolos, nivel_bits);

for i = 1:num_simbolos
    indice_inicial = (i - 1) * nivel_bits + 1;
    indice_final = i * nivel_bits;
    simbolo = senal_normalizada(indice_inicial:indice_final);
    bits_cuantificados(i, :) = simbolo >= 0.5; % Cuantificación binaria (0 o 1)
end

b = reshape(bits_cuantificados', [], 1); % Convertir matriz en vector columna

bn=num2str(b);
bx=['La cadena binaria aleatoria es ',bn']

%% SEÑALES SENO Y COSENO
in=cos(f*t); % Coseno
qd=sin(f*t); % Seno
cp=[];sp=[];
A=1/sqrt(2);         % Amplitud seno y coseno
mod1=[];mod2=[];bit=[];
m=A*ones(1,sam);

%% MODULACIÓN QPSK
for j=1:2:length(b)
    eb=bin2int(b(j:j+1));
    switch eb 
        case 0
          m1=-m;
          m2=m1;
          se=zeros(1,sam);
     case 1
          m1=m;
          m2=-m;
          se=[zeros(1,sam/2),ones(1,sam/2)];
     case 2
          m1=-m;
          m2=m;
          se=[ones(1,sam/2),zeros(1,sam/2)];
     case 3
          m1=m;
          m2=m;
          se=ones(1,sam);
        otherwise 
    end 
    cp=[cp,m1]; %#ok<AGROW> 
    sp=[sp,m2]; %#ok<AGROW> 
    mod1=[mod1,in]; %#ok<AGROW> 
    mod2=[mod2,qd]; %#ok<AGROW> 
    bit=[bit,se]; %#ok<AGROW> 
end



%% RESULTADOS MODULACIÓN
qpsk_complex=complex(cp,sp);
qpsk=cp.*mod1+sp.*mod2;
figure(1); % SEÑAL CODIFICADA (SEÑAL DIGITAL)
subplot(211)
plot(bit,'Color','[0.6350 0.0780 0.1840]','LineWidth',1.5');
grid on;
title('Random Binary Signal');
xlabel('Tiempo');ylabel('Amplitud');

subplot(212) %SEÑAL MODULADA TX
plot(qpsk,'k','LineWidth',1.5');
grid on;
title('QPSK Output Waveform');
xlabel('Tiempo');ylabel('Amplitud');

%% ENTROPIA DE LA FUENTE 


probabilidades = zeros(1, 2^nivel_bits); % Inicializar vector de probabilidades

for i = 1:num_simbolos
    simbolo = bi2de(bits_cuantificados(i, :)) + 1; % Convertir binario a decimal y ajustar índice
    probabilidades(simbolo) = probabilidades(simbolo) + 1; % Incrementar contador para el símbolo correspondiente
end

probabilidades = probabilidades / num_simbolos; % Calcular probabilidades
probabilidades = probabilidades(probabilidades > 0); % Eliminar ceros en caso de existir

H_transmisor = -sum(probabilidades .* log2(probabilidades)); % Calcular entropía

disp("Entropía en el transmisor: " + num2str(H_transmisor) + " bits");

%% Add Noise 
%Ploteo Modulada y con Ruido
snr=10;
rxSig = awgn(qpsk,snr);
figure(2)
subplot(211),plot(qpsk,'b');grid on;
xlim([0 6700])
title('QPSK Output Waveform');
xlabel('Tiempo');ylabel('Amplitud');
subplot(212),plot(rxSig,'g');grid on;
xlim([0 6700])
title('QPSK Output with Noise Waveform');
xlabel('Tiempo');ylabel('Amplitud');
%% Constellation and Decision Diagram
cp_snr=awgn(cp,snr);
sp_snr=awgn(sp,snr);
qpsk_snr=cp_snr.*mod1+sp_snr.*mod2;
z=complex(cp,sp);
Ireal=cp(100:100:length(cp));
Qimag=sp(100:100:length(cp));
zplot= complex(Ireal,Qimag);
hold off
Ir_awgn=cp_snr(100:100:length(cp_snr));
Qi_awgn=sp_snr(100:100:length(cp_snr));
zplot_awgn=complex(Ir_awgn,Qi_awgn);
%Region decision - Diagrama Constelacions antes del SNR
figure(3)
scatter(Ireal,Qimag,'m','o');
xlabel('Phase')
ylabel('Quadrature')
title('Diagrama de constelación antes del canal')
grid on
hold on
au=[-1.5 1.5];
bu=au-au;
plot(au,bu,'k');
plot(bu,au,'k');
hold off
%Diagrama Constelacion Despues SNR
figure(4)
scatter(Ir_awgn,Qi_awgn,'m','o');
xlabel('Phase')
ylabel('Quadrature')
title('Diagrama de constelación despues de AWGN')
grid on
hold on
au=[-1.5 1.5];
bu=au-au;
plot(au,bu,'k');
plot(bu,au,'k');
hold off

%Sobreposicion de los diagramas
figure(5)
scatter(Ireal,Qimag,'r','o');
hold on
scatter(Ir_awgn,Qi_awgn,'g','o');
xlabel('Phase')
ylabel('Quadrature')
title('Diagrama de constelación')
grid on
hold on
au=[-1.5 1.5];
bu=au-au;
plot(au,bu,'k');
plot(bu,au,'k');
hold off

%% CALCULO BER PREDETECCION
snrdB=0:1:snr;
snrlin=10.^(snrdB/10);
BERteorico = 0.5*erfc(sqrt(snrlin));
figure(6)
semilogy(snrdB,BERteorico);
hold on;
for i = 1:length(snrdB)
    Qpsk_BER=awgn(z,i);
    in_Rx=real(Qpsk_BER);
    qd_Rx=imag(Qpsk_BER);
    
    cp_Tx=(cp > 0);
    sp_Tx=(sp > 0);
    
    cp_Rx=(in_Rx > 0);
    sp_Rx=(qd_Rx > 0);
    
    in_error=(cp_Rx~=cp_Tx);
    qd_error=(sp_Rx~=sp_Tx);

    bit_error=sum(in_error)+sum(qd_error); %
    BER(i)=sum(bit_error)/(2*length(Qpsk_BER)); 
   
       
end
semilogy(snrdB,BER,'o')
grid on
legend('BER Teórico','BER Simulado')
title('BER simulado vs BER teórico')
hold off
%% Entropia en el receptor


H_receptor=entropy(rxSig);
disp("Entropía en el receptor: " + num2str(H_receptor) + " bits");


%% Demodulacion 
for i=1:length(b)/2
    %%Detector de Fase coherente
    Z_in=qpsk((i-1)*length(t)+1:i*length(t)).*cos(f*t); 
    Z_in_intg=(trapz(t,Z_in))*(2*f);
    if(Z_in_intg>0) 
        r_in=1;
    else
        r_in=0; 
    end
%Detector de cuadratura coherente
    Z_qd=qpsk((i-1)*length(t)+1:i*length(t)).*sin(f*t);
    Z_qd_intg=(trapz(t,Z_qd))*(2*f);
    if (Z_qd_intg>0)
        r_qd=1;
    else
        r_qd=0; 
    end    
end
rxData=[];
for i=1:length(b)/2
    %Detector coherente en fase
    Z_in=rxSig((i-1)*length(t)+1:i*length(t)).*cos(f*t); 
    Z_in_intg=(trapz(t,Z_in))*(2*f);
    if(Z_in_intg>0) 
        r_in=1;
    else
        r_in=0; 
    end
%Detector de Cuadratura Coherente
    Z_qd=rxSig((i-1)*length(t)+1:i*length(t)).*sin(f*t);
    Z_qd_intg=(trapz(t,Z_qd))*(2*f);
    if (Z_qd_intg>0)
        r_qd=1;
    else
        r_qd=0; 
    end
    rxData=[rxData  r_qd  r_in];
end

Rx=[];
for i=1:length(b)
    Rx=[Rx rxData(i)*ones(1,sam)];
end
figure(7)
subplot(211), plot(bit,'k','LineWidth',1'),ylim([-0.2 1.2])
title('Bits de información transmitidos')
subplot(212), plot(Rx, 'r','LineWidth',1'),ylim([-0.2 1.2])
title('Bits de información recibidos')



