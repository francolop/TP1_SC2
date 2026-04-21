% =========================================================
% Alumno  : Lopez Franco
% Profesor: Pucheta Julian
% Materia : Control 2
% =========================================================
% Item [4]: Motor de DC por metodo de euler
% =========================================================

clear; clc; close all;

%% 1) PARÁMETROS DEL SISTEMA
La = 366e-6;     % Inductancia de armadura [H]
Jm = 5e-9;       % Inercia del rotor [kg*m^2]
Ra = 55.6;       % Resistencia [Ohm]
Bv = 0;          % Fricción viscosa
Kt = 6.49e-3;    % Constante de torque [Nm/A]
Ke = 6.53e-3;    % Constante FEM [V/(rad/s)]

%% 2) MODELO EN ESPACIO DE ESTADOS

M = [ -Ra/La    -Ke/La    0;
       Kt/Jm    -Bv/Jm    0;
       0         1        0 ];

N = [ 1/La   0;
      0    -1/Jm;
      0     0 ];

%% 3) PARÁMETROS DE SIMULACIÓN
h  = 10e-6;       % paso temporal
Tf = 5;          % tiempo final

time = 0:h:Tf;
Ns = length(time);

%% 4) VARIABLES

states = zeros(3, Ns);   % [ia; wr; theta]
inputs = zeros(2, Ns);   % [Va; TL]

% Señales de entrada
Va = 12;
TL = 0;

inputs(1,:) = Va;
inputs(2,:) = TL;

% Condición inicial
states(:,1) = [0; 0; 0];

%% 5) INTEGRACIÓN (EULER)

for k = 1:Ns-1
    dx = M*states(:,k) + N*inputs(:,k);
    states(:,k+1) = states(:,k) + h*dx;
end

%% 6) VARIABLES DE SALIDA

corriente = states(1,:);
velocidad = states(2,:);
angulo    = states(3,:);

%% 7) GRÁFICOS

figure

subplot(3,1,1)
plot(time, corriente)
title('Corriente de armadura')
ylabel('A')
grid on

subplot(3,1,2)
plot(time, velocidad)
title('Velocidad angular')
ylabel('rad/s')
grid on

subplot(3,1,3)
plot(time, angulo)
title('Posición angular')
ylabel('rad')
xlabel('Tiempo [s]')
grid on