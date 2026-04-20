% =========================================================
% Alumno  : Lopez Franco
% Profesor: Pucheta Julian
% Materia : Control 2
% =========================================================
% Item [1]:  Obtener simulaciones que permitan estudiar la dinámica del sistema
% =========================================================


clear all;
% Parámetros de diseño
R = 2200; L = 10e-6; Cap = 500e-6; vin = 12;   
% Punto de operación
Il(1) = 0;
Vcl(1)= 0;
y(1)  = 0;
Xop   = [0 0]';          
x     = [Il(1) Vcl(1)]'; % variables de estado

% Espacio de estados
A=[-R/L -1/L; 1/Cap 0];   % matriz de estados
B=[1/L; 0];               % matriz de entrada
C=[R 0];                  % matriz de salida
D=[0];                    % matriz de transmision directa

% Cálculo de la funcion de transferencia
[numG, denG] = ss2tf(A,B,C,D);
% Función de transferencia del sistema
G = tf(numG, denG)
% Polos de la FT
poles = roots(denG)


tR   = log(0.95)/poles(1)    % dinámica rápida

tint = tR/10;                 % tiempo de integración (10 veces menor)

tsim = 2e-3                   
step = tsim/tint;
t    = linspace(0, tsim, step);
u    = linspace(0, 0, step);    % vector inicial de entrada (todos 0s)

%%
% flag para controlar switcheo de la entrada 
swTime = 0;     
z      = size(t);
z      = z(2)/2;

for i=1:step-1
 swTime = swTime+1;
 if (swTime >= z)
     swTime = 0;
     vin = vin*(-1);
 end
 
 u(i) = vin;
 
 %Variables de sistema lineal 
 xp = A*(x-Xop)+B*u(i); % calcula derivadas de las VE
 x = x + xp*tint;       % actualiza el valor de las VE para el próximo
                        % paso por integración por Euler
 Y = C*x;               % calcula salida
 % Actualización de las variables de interés para los nuevos instantes de
 % tiempo
 y(i+1)=Y(1);           % actualiza valor de la salida
 Il(i+1)=x(1);          % actualiza valor de la corriente en el circuito
 Vcl(i+1)=x(2);         % actualiza valor de la tensión en C
end

%%
% Gráficas

% Gráfico de Il, Vc y Vin
figure(1)
subplot(4,1,1);
plot(t,Il, 'b' );title('Corriente , i'); grid on; 
subplot(4,1,2);
plot(t,Vcl, 'r' );title('Tensión Capacitor , v_c');grid on
subplot(4,1,3); 
plot(t,u, 'm' );title('Tensión de Entrada, u');grid on
subplot(4,1,4);
plot(t,y, 'm' );title('Tensión de Salida, v_r');grid on
