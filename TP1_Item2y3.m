% =========================================================
% Alumno  : Lopez Franco
% Profesor: Pucheta Julian
% Materia : Control 2
% =========================================================
% Item [2]: Identificacion de R, L y C mediante respuesta al escalon
% Item [3]: Validacion superponiendo curvas de corriente desde t=0.05s
% =========================================================

clear; close all; clc;
%  IMPORTACION DE DATOS
[time_v, i_v, Vc_v, Ve_v] = import_data("Curvas_Medidas_RLC_2026.xls", 1, 2, 2001);

Ts = time_v(2) - time_v(1);  

%  ITEM [2] — METODO DE LA RESPUESTA AL ESCALON

% --- Deteccion automatica del retardo (delay) ---
idx_delay = find(Ve_v > max(Ve_v)/2, 1, 'first');
delay     = time_v(idx_delay);

% --- Amplitud del escalon ---
StepAmplitude = max(Ve_v);

% --- Tiempos de muestreo dentro del transitorio ---
t1 = 0.005;
t2 = 2 * t1;
t3 = 3 * t1;

% Busqueda de indices por minima distancia (robusto ante cualquier offset)
[~, idx1] = min(abs(time_v - (delay + t1)));
[~, idx2] = min(abs(time_v - (delay + t2)));
[~, idx3] = min(abs(time_v - (delay + t3)));

y1 = Vc_v(idx1);
y2 = Vc_v(idx2);
y3 = Vc_v(idx3);

% --- Normalizacion ---
K  = 1;
k1 = (1/StepAmplitude) * y1/K - 1;
k2 = (1/StepAmplitude) * y2/K - 1;
k3 = (1/StepAmplitude) * y3/K - 1;

% --- Discriminante ---
be = 4*k1^3*k3 - 3*k1^2*k2^2 - 4*k2^3 + k3^2 + 6*k1*k2*k3;

if be >= 0
    alfa1 = (k1*k2 + k3 - sqrt(be)) / (2*(k1^2 + k2));
    alfa2 = (k1*k2 + k3 + sqrt(be)) / (2*(k1^2 + k2));
else
    warning('be < 0 (sistema subamortiguado): se usa alfa1 = alfa2.');
    alfa_r = (k1*k2 + k3) / (2*(k1^2 + k2));
    alfa1  = alfa_r;
    alfa2  = alfa_r;
end

% --- Constantes de tiempo ---
T1 = -t1 / log(alfa1);
T2 = -t1 / log(alfa2);

% --- Funcion de transferencia identificada ---
G_id = tf(1, conv([T1 1], [T2 1]));

fprintf('[Item 2] Funcion de transferencia identificada:\n');
G_id   % muestra la TF en formato compacto

% --- Identificacion de R, L y C ---
i_post = i_v(idx_delay:end);       
R      = StepAmplitude / max(i_post); 
C_val = (T1 + T2) / R;              % C*R = T1+T2
L_val = (T1 * T2) / C_val;          % C*L = T1*T2

fprintf('[Item 2] Parametros identificados del circuito RLC:\n');
fprintf('  R = %.4f  Ohm\n', R);
fprintf('  C = %.4e  F\n',   C_val);
fprintf('  L = %.4e  H\n\n', L_val);

% --- Verificacion con funcion de transferencia teorica (simbolica) ---
syms Rs Ls Cs s_sym
A_sym   = [-Rs/Ls, -1/Ls; 1/Cs, 0];
B_sym   = [1/Ls; 0];
C_sym   = [0 1];
G_sym   = simplify((C_sym * adjoint(s_sym*eye(2) - A_sym) * B_sym) / ...
                    det(s_sym*eye(2) - A_sym));
fprintf('[Verificacion] G teorica (simbolica):\n  ');
disp(G_sym);





%  ITEM [3] — VALIDACION POR SIMULACION (Euler hacia adelante)
A_m = [-R/L_val,  -1/L_val;
        1/C_val,   0      ];
B_m = [1/L_val; 0];

x   = [0; 0];   % condiciones iniciales nulas
N   = length(time_v);
Vc_sim = zeros(1, N);
I_sim  = zeros(1, N);

for k = 1:N
    I_sim(k)  = x(1);
    Vc_sim(k) = x(2);
    x = x + (A_m*x + B_m*Ve_v(k)) * Ts;
end

% Indice desde t = 0.05 s para comparacion de corriente (Item 3)
[~, idx_05] = min(abs(time_v - 0.05));



%  GRAFICAS
fz  = 14;
lw1 = 1.8;
lw2 = 1.5;

% --- Figura 1: Tension en el capacitor ---
figure('Name', 'Item 2 - Tension en el Capacitor');
plot(time_v, Vc_sim, 'b',   'LineWidth', lw1); hold on;
plot(time_v, Vc_v,   'r--', 'LineWidth', lw2);
plot(time_v, Ve_v,   'k:',  'LineWidth', 1.2);
xline(delay, 'g--', 'Inicio escalon', 'LabelVerticalAlignment', 'bottom');
grid on;
% Titulo sin latex para evitar conflicto con corchetes
title('Item 2 - Tension en el Capacitor', 'FontSize', fz+1);
ylabel('Vc [V]', 'Interpreter', 'latex', 'FontSize', fz);
xlabel('t [s]',  'Interpreter', 'latex', 'FontSize', fz);
legend('Simulada (R,L,C identificados)', 'Medida', 'Entrada Ve', ...
       'Interpreter', 'latex', 'FontSize', fz-2, 'Location', 'southeast');

% --- Figura 2: Corriente desde t = 0.05 s ---
figure('Name', 'Item 3 - Corriente desde t=0.05s');
plot(time_v(idx_05:end), I_sim(idx_05:end), 'b',   'LineWidth', lw1); hold on;
plot(time_v(idx_05:end), i_v(idx_05:end),   'r--', 'LineWidth', lw2);
grid on;
title('Item 3 - Validacion: Corriente i(t) desde t = 0.05 s', 'FontSize', fz+1);
ylabel('i [A]', 'Interpreter', 'latex', 'FontSize', fz);
xlabel('t [s]', 'Interpreter', 'latex', 'FontSize', fz);
xlim([0.05, time_v(end)]);   % mostrar todo el rango disponible desde 0.05s
legend('Simulada (R,L,C identificados)', 'Medida', ...
       'Interpreter', 'latex', 'FontSize', fz-2, 'Location', 'northeast');
   
   
   


function [time_v, i_t, Vc_t, Ve_t] = import_data(workbookFile, sheetName, startRow, endRow)
%IMPORT_DATA  Lee las 4 columnas de datos del archivo Excel.

    if nargin < 2 || isempty(sheetName), sheetName = 1; end
    if nargin < 4
        startRow = 2;
        endRow   = 2001;
    end

    opts = spreadsheetImportOptions("NumVariables", 4);
    opts.Sheet         = sheetName;
    opts.DataRange     = "A" + startRow(1) + ":D" + endRow(1);
    opts.VariableNames = ["time_v", "i_t", "Vc_t", "Ve_t"];
    opts.VariableTypes = ["double", "double", "double", "double"];

    tbl = readtable(workbookFile, opts, "UseExcel", false);

    for k = 2:length(startRow)
        opts.DataRange = "A" + startRow(k) + ":D" + endRow(k);
        tbl = [tbl; readtable(workbookFile, opts, "UseExcel", false)]; %#ok<AGROW>
    end

    time_v = tbl.time_v;
    i_t    = tbl.i_t;
    Vc_t   = tbl.Vc_t;
    Ve_t   = tbl.Ve_t;
end
