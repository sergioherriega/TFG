function [Rx]=correfMcss(M,x,L,W)
% correfMcss realiza la correlacion eficiente con un conjunto de M secuencias complementarias
%
% [Rx]=correfMcss(M,x,L,W);
% M:Numero de elementos del conjunto.
% x:Secuencia de entrada.
% L:longitud de las secuencias complementarias (M^N).
% W: Numero natural menor que L cuyo equivalente binario de n*N bits representa 
% a las semillas que se utilizaron en la generacion de las secuencias complementarias
% Rx: Salidas de la señal correlada con cada secuencia conjunto generada
% con la semilla W.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% === CAMBIO MÍNIMO PARA MATLAB CODER ===
% Asegura que Rx esté definido incluso si se entra en algún "return" de validación
Rx = zeros(M, length(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chequeo de los argumentos de entrada %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if M<2
    disp('Error: M debe ser mayor o igual que 2');
    return
end
if mod(log(M)/log(2),1)~=0
    disp('Error: M debe ser potencia de 2');
    return
end
N=log(L)/log(M);
if mod(N,1)~=0
    disp('Error: L debe ser igual a una potencia de M');
    return
end
if W>=L
    disp('Error: W debe ser menor que L');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Definicion de variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=log2(M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conversión de W a vector binario ANTES de indexarlo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnr = dec2bin(W,n*N);
W   = (bnr=='1')-(bnr~='1'); % transforma los ceros en -1s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generacion de coeficientes para 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:N
    coef=[1 W(n*j-n+1);1 -W(n*j-n+1)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generacion de coeficientes para cada etapa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coefn=zeros(M,M,N);
coef=[1 W(1);1 -W(1)];
coefn(:,:,1)=coef;

for j=2:N
    coef=[1 W(n*j-n+1);1 -W(n*j-n+1)];
    coefn(:,:,j)=coef;
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inicializacion de variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Valores acumulados (Buffers de memoria)
R = zeros(M-1, L);               % buffer fijo
% Inicializacion de las salidas
Rx = zeros(M, length(x));        % salida fija: M x length(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Algoritmo principal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Itera la señal de entrada 
for i=1:length(x) 
    ind = zeros(1, M-1);
    Rdt = zeros(M, N);
    
  %Etapa 1
    retardo = L/M;   % primera etapa de retardos
    for j=1:M-1
        ind(j) = 1 + retardo*(M-j);
    end

    % Actualiza retardos (shift register)
    for k=1:M-1
        R(k,2:end) = R(k,1:end-1);
        R(k,1)     = x(i);
    end

    % Correlación jerárquica
    Rdt(M,1) = coefn(M,M,1)*R(M-1,1);
    for k=1:M-1
        Rdt(M,1) = Rdt(M,1) + coefn(M,k,1)*R(k,ind(k));
    end

    for k=1:M-1
        Rdt(k,:) = coefn(k,M,1)*R(M-1,1);
        for m=1:M-1
            Rdt(k,:) = Rdt(k,:) + coefn(k,m,1)*R(m,ind(m));
        end
    end
    
    % ACTUALIZACIÓN DE RETARDOS (ESCALAR)
    for k=1:M-1
        R(k,ind(k)) = Rdt(k,1);
    end
    
  %Etapas 2 hasta N
    for j=2:N
        retardo = retardo/M;
        indn = zeros(1, M-1);
        for k=1:M-1
            indn(k) = ind(k) + retardo*(M-k);
        end
        ind = indn;

        Rdt(M,j) = coefn(M,M,j)*Rdt(M,j-1);
        for k=1:M-1
            Rdt(M,j) = Rdt(M,j) + coefn(M,k,j)*R(k,ind(k));
        end

        for k=1:M-1
            Rdt(k,:) = coefn(k,M,j)*Rdt(M,j-1);
            for m=1:M-1
                Rdt(k,:) = Rdt(k,:) + coefn(k,m,j)*R(m,ind(m));
            end
        end
        
        for k=1:M-1
            R(k,ind(k)) = Rdt(k,j);
        end
    end
    
  % Resultado final
    for j=1:M
        if j==M
            Rx(j,i) = Rdt(M,N);
        else
            Rx(j,i) = R(j,ind(j));
        end
    end
end
end