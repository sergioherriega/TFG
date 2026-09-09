function [Ry_gk] = CorrLSmiracle_VELSC(g, seed, N, Lo, K, Code)
% CorrLS_modificado: Implementa el correlador adaptado al segundo diagrama
% Entradas:
%   g: Señal de entrada
%   seed: Semilla para secuencia complementaria
%   N: Número de secuencias CCS
%   Lo: Longitud de secuencias CCS
%   K: Número de códigos
%   Code: Índice de fila en matriz de Hadamard

% Parámetros iniciales
M = 2; % Golay
W0 = Lo - 1; % Longitud de zona libre de interferencia
p = K / 2; % Mitad de los códigos
pii = [1, 1, 2, 2]; % Permutación
H = [1, 1, 1, 1; 1, -1, 1, -1; 1, 1, -1, -1; 1, -1, -1, 1]; %Matriz de Hadamard 

% Obtener códigos complementarios desde las semillas
if seed < Lo/M
    seed_1 = Lo/M + seed;
else
    seed_1 = seed - Lo/M;
end

% BLOQUE 1 2: correlación con códigos complementarios
[setL0] = correfMcss(N, [g zeros(1,Lo)], Lo, seed);
[setL1] = correfMcss(N, [g zeros(1,Lo)], Lo, seed_1);

% Separación de componentes C y S
Lrecep = numel(g) + Lo;   % <-- CAMBIO MÍNIMO (antes era 2*Lo)
C = zeros(2, Lrecep);
S = zeros(2, Lrecep);
C(1,:) = setL0(1,:); C(2,:) = setL1(1,:);
S(1,:) = setL0(2,:); S(2,:) = setL1(2,:);
% Lrecep ya es fijo arriba

% Preparar salida: tamaño máximo posible considerando retardos
maxDelay = max((K/2 + (4 - 1)) * Lo + W0, (K/2 - 1) * Lo);
outLen = maxDelay + Lrecep;
Ry_gk  = zeros(1, outLen);

% BLOQUES 3 4: retardos, ponderaciones, combinación cruzada
for i = 1:p
    % Retardos 
    d1 = (K/2 + (4 - i)) * Lo + W0;
    d2 = (K/2 - i) * Lo;

    % Inicializar 
    x0 = zeros(1, outLen);
    y0 = zeros(1, outLen);

    % Aplicar ponderaciones Hadamard
    x0(d1+1 : d1+Lrecep) = H(Code,i) * C(pii(i), :);
    y0(d2+1 : d2+Lrecep) = H(Code,i) * S(pii(i), :);

    % Ventanas 
    WT = 1;
    W1 = 1;

    % Aplicación final
    xT = x0;
    yT = WT * y0;

    Ry_gk = Ry_gk + xT + W1 * yT;
end
% === Graficar señales internas ===
% figure(2); plot(x0); grid on; title('Señal x0'); xlabel('Desplazamiento (k)'); ylabel('Valor de Correlación');
% figure(3); plot(y0); grid on; title('Señal y0'); xlabel('Desplazamiento (k)'); ylabel('Valor de Correlación');
% figure(4); plot(xT); grid on; title('Señal xT'); xlabel('Desplazamiento (k)'); ylabel('Valor de Correlación');
% figure(5); plot(yT); grid on; title('Señal yT'); xlabel('Desplazamiento (k)'); ylabel('Valor de Correlación');

end