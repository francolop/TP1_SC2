% =========================================================
% Alumno  : Lopez Franco
% Profesor: Pucheta Julian
% Materia : Control 2
% =========================================================
% Item [6]: PID en tiempo discreto 
% =========================================================


clc; clear all; close all;

%% PARÁMETROS DEL MOTOR 
La  = 7.0274e-4;   % [H]
Ra  = 28.13;       % [Ohm]
Km  = 0.0605;      % [V·s/rad]
Ki_motor = 0.01162;% [N·m/A]
Jm  = 2.0628e-9;   % [kg·m²]
Bm  = 0;           % [N·m·s/rad]

%% ESPACIO DE ESTADOS CONTINUO
% x = [ia; wr; theta]  |  u = [Va; TL]
A = [-Ra/La,        -Km/La,  0;
      Ki_motor/Jm,  -Bm/Jm,  0;
      0,             1,       0];

Bmat = [1/La,   0;
        0,     -1/Jm;
        0,      0];

C = eye(3);
D = zeros(3,2);

sys_c = ss(A, Bmat, C, D);

%% DISCRETIZACIÓN
Ts    = 1e-4;           
sys_d = c2d(sys_c, Ts, 'zoh');
Ad    = sys_d.A;
Bd    = sys_d.B;

%% TIEMPO DE SIMULACIÓN
tf = 40;
t  = 0:Ts:tf;
N  = length(t);

%% REFERENCIA Y PERTURBACIÓN
theta_ref = ones(1,N);   

% Torque escalado a los parámetros del motor
% Torque nominal: T = Ki * (Va/Ra)
TL_amp = Ki_motor * (12/Ra) * 0.5;  
fprintf('TL_amp = %.6e N·m\n', TL_amp);

TL = zeros(1,N);
TL(t >= 19 & t < 27) = TL_amp;

%% PID DISCRETO INCREMENTAL
% u(k) = u(k-1) + A0*e(k) + B0*e(k-1) + C0*e(k-2)
Kp = 10;
Ki = 200;
Kd = 5e-5;

A0 = (2*Kp + Ki*Ts + 2*Kd/Ts) / 2;
B0 = (-2*Kp + Ki*Ts - 4*Kd/Ts) / 2;
C0 = Kd / Ts;

fprintf('\n===== PID DISCRETO =====\n');
fprintf('Kp = %.4f\n', Kp);
fprintf('Ki = %.4f\n', Ki);
fprintf('Kd = %.2e\n', Kd);
fprintf('A0 = %.4e\n', A0);
fprintf('B0 = %.4e\n', B0);
fprintf('C0 = %.4e\n', C0);

u_max =  12;
u_min = -12;

%% INICIALIZACIÓN
x     = zeros(3,N);
u     = zeros(1,N);
e     = zeros(1,N);

x(:,1) = [0;0;0];
e(1)   = theta_ref(1) - x(3,1);
e(2)   = theta_ref(2) - x(3,2);

%% SIMULACIÓN EN LAZO CERRADO
for k = 3:N-1
    e(k) = theta_ref(k) - x(3,k);
    
    % PID incremental
    u(k) = u(k-1) + A0*e(k) + B0*e(k-1) + C0*e(k-2);
    
    % Saturación
    u(k) = min(max(u(k), u_min), u_max);
    
    % Actualización de estados
    x(:,k+1) = Ad*x(:,k) + Bd*[u(k); TL(k)];
end
e(N) = theta_ref(N) - x(3,N);

ia    = x(1,:);
wr    = x(2,:);
theta = x(3,:);

%% RESULTADOS EN CONSOLA
fprintf('\n===== RESULTADOS =====\n');
fprintf('Error final de ángulo  = %.6e rad\n', theta_ref(end)-theta(end));
fprintf('Máx corriente          = %.6e A\n',   max(abs(ia)));
fprintf('Máx velocidad          = %.6e rad/s\n', max(abs(wr)));
fprintf('Máx tensión control    = %.6e V\n',   max(abs(u)));
fprintf('Máx ángulo             = %.6e rad\n', max(theta));
fprintf('Mín ángulo             = %.6e rad\n', min(theta));

%% FIGURA 1 - Vista completa
figure(1);
subplot(4,1,1);
plot(t, theta_ref, 'k--', 'LineWidth', 1.2); hold on;
plot(t, theta, 'b', 'LineWidth', 1.3);
grid on; ylabel('\theta [rad]');
title('Seguimiento del ángulo');
legend('\theta_{ref}','\theta','Location','best');
ylim([0.8 1.2]);

subplot(4,1,2);
plot(t, wr, 'r', 'LineWidth', 1.3);
grid on; ylabel('\omega_r [rad/s]');
title('Velocidad angular');

subplot(4,1,3);
plot(t, ia, 'm', 'LineWidth', 1.3);
grid on; ylabel('i_a [A]');
title('Corriente de armadura');

subplot(4,1,4);
plot(t, u,  'b', 'LineWidth', 1.3); hold on;
plot(t, TL, 'k--', 'LineWidth', 1.0);
grid on; xlabel('Tiempo [s]');
ylabel('u [V] / T_L');
title('Acción de control y perturbación');
legend('u[k]','T_L','Location','best');
sgtitle('PID discreto - Parámetros inciso 5');

%% FIGURA 2 - Zoom transitorio inicial
figure(2);
idx_ini = t <= 0.5;

subplot(4,1,1);
plot(t(idx_ini), theta_ref(idx_ini), 'k--', 'LineWidth', 1.2); hold on;
plot(t(idx_ini), theta(idx_ini), 'b', 'LineWidth', 1.3);
grid on; ylabel('\theta [rad]');
title('Zoom inicial - Ángulo');
legend('\theta_{ref}','\theta');

subplot(4,1,2);
plot(t(idx_ini), wr(idx_ini), 'r', 'LineWidth', 1.3);
grid on; ylabel('\omega_r [rad/s]');
title('Zoom inicial - Velocidad');

subplot(4,1,3);
plot(t(idx_ini), ia(idx_ini), 'm', 'LineWidth', 1.3);
grid on; ylabel('i_a [A]');
title('Zoom inicial - Corriente');

subplot(4,1,4);
plot(t(idx_ini), u(idx_ini), 'b', 'LineWidth', 1.3);
grid on; xlabel('Tiempo [s]'); ylabel('u [V]');
title('Zoom inicial - Acción de control');
sgtitle('Zoom transitorio inicial');

%% FIGURA 3 - Zoom perturbación
figure(3);
idx_pert = (t >= 17.5) & (t <= 27.5);

subplot(4,1,1);
plot(t(idx_pert), theta_ref(idx_pert), 'k--', 'LineWidth', 1.2); hold on;
plot(t(idx_pert), theta(idx_pert), 'b', 'LineWidth', 1.3);
grid on; ylabel('\theta [rad]');
title('Zoom perturbación - Ángulo');
legend('\theta_{ref}','\theta');

subplot(4,1,2);
plot(t(idx_pert), wr(idx_pert), 'r', 'LineWidth', 1.3);
grid on; ylabel('\omega_r [rad/s]');
title('Zoom perturbación - Velocidad');

subplot(4,1,3);
plot(t(idx_pert), ia(idx_pert), 'm', 'LineWidth', 1.3);
grid on; ylabel('i_a [A]');
title('Zoom perturbación - Corriente');

subplot(4,1,4);
plot(t(idx_pert), u(idx_pert),  'b', 'LineWidth', 1.3); hold on;
plot(t(idx_pert), TL(idx_pert), 'k--', 'LineWidth', 1.0);
grid on; xlabel('Tiempo [s]'); ylabel('u [V] / T_L');
title('Zoom perturbación - Control');
legend('u[k]','T_L');
sgtitle('Respuesta frente a la perturbación');