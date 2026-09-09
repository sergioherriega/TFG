%Correlador de secuencias LS, usando funcion CorreladorLSmiracle_VELSC
%8HC

clc;
close all;
clear;


%% README

%% Cargar Secuencia a correlar
load("Secuencia_LS_1151.mat");

%% Variables
Code = 3; %
g = G(Code,:); % Selecciona secuencia 1 de grupo de secuencias LS generadas.
seed = 0; % Semilla 0, se define en la generación.
N = 2; % Número de secuencias a la salida del correlador (Normal y Ortogonal)
W = 127; % ventana libre de interferncia
K = 8; % Número de codigos de una familia
Lo = (L - W)/K; % Longitud de la secuencia CCS

%% Funcion Generadora
[Crgk] = CorrLSmiracle_VELSC(g,seed,N,Lo,K,Code);

%% Visualización --> en este caso no las queremos porque el MatlabCoder no lo tolera asiqeu comento TODAS las representaciiones
Crgk_col = Crgk(1,:)';
figure(1)
hold on
grid on
plot(Crgk(1,:));
grid on; 
title('Correlación Cruzada de la Secuencia LS usando Correlador VELSC');
xlabel('Desplazamiento (k)');
ylabel('Valor de Correlación');
