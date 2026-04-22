% =========================================================
% Alumno  : Lopez Franco
% Profesor: Pucheta Julian
% Materia : Control 2
% =========================================================
% Item [5]: Otencion del modelo del motor mediante curvas
% =========================================================





clc; clear; close all;

% load data
raw = readmatrix('Curvas_Medidas_Motor_2026');
time  = raw(:,1);
speed = raw(:,2);
curr  = raw(:,3);
volt  = raw(:,4);
load_ = raw(:,5);
ts    = mean(diff(time));

% smooth signals
spd_s = smoothdata(speed, 'sgolay', 21);
cur_s = smoothdata(curr,  'sgolay', 21);

% find step
dv      = diff(volt);
stepIdx = find(dv > 0.5*max(dv), 1, 'first') + 1;
if isempty(stepIdx)
    error('Step not found.');
end
t0  = time(stepIdx);
Ei  = max(volt);
fprintf('Step at t0 = %.4f s,  Ei = %.4f V\n', t0, Ei);

% identification window
t_a   = t0;
t_b   = min(t0 + 8, time(end));
mask  = (time >= t_a) & (time <= t_b);
t_w   = time(mask);
spd_w = spd_s(mask);
cur_w = cur_s(mask);
v_w   = volt(mask);

% initial / final values
W0 = mean(spd_s(max(1, stepIdx-20) : stepIdx-1));
I0 = mean(cur_s(max(1, stepIdx-20) : stepIdx-1));
eI = find(mask, 1, 'last');
Wf = mean(spd_s(max(1, eI-20) : eI));
If = mean(cur_s(max(1, eI-20) : eI));

Kw = (Wf - W0) / Ei;
Ki = (If - I0) / Ei;
fprintf('Kw = %.8f  Ki = %.8f\n', Kw, Ki);

% Chen points search
fprintf('Searching Chen points - speed...\n');
resW = findChenPoints(time, spd_s, W0, Ei, Kw, t0, t_a, t_b, volt, ...
    t0+0.6, t0+4.0, 0.05, 0.15);

fprintf('Searching Chen points - current...\n');
resI = findChenPoints(time, cur_s, I0, Ei, Ki, t0, t_a, t_b, volt, ...
    t0+0.4, t0+3.5, 0.05, 0.15);

% print results
fprintf('\n--- SPEED TF ---\n');
fprintf('Chen pts: [%.4f %.4f %.4f] s\n', resW.tpts_abs);
fprintf('T1=%.6f  T2=%.6f  T3=%.6f  NRMSE=%.6f\n', resW.T1, resW.T2, resW.T3, resW.nrmse);

fprintf('\n--- CURRENT TF ---\n');
fprintf('Chen pts: [%.4f %.4f %.4f] s\n', resI.tpts_abs);
fprintf('T1=%.6f  T2=%.6f  T3=%.6f  NRMSE=%.6f\n', resI.T1, resI.T2, resI.T3, resI.nrmse);

Gw = resW.G;
Gi = resI.G;

% simulate TF
t_sim = time - time(1);
spd_tf = lsim(Gw, volt, t_sim) + W0;
cur_tf = lsim(Gi, volt, t_sim) + I0;

% extract motor parameters
[nW, dW] = tfdata(Gw, 'v');
[nI, dI] = tfdata(Gi, 'v');
nW = nW/dW(end); dW = dW/dW(end);
nI = nI/dI(end); dI = dI/dI(end);
tol = 1e-12;
nW = nW(find(abs(nW)>tol,1):end);
nI = nI(find(abs(nI)>tol,1):end);

Km_e = nW(2);
J_e  = nI(1);
B_e  = nI(2);
L_e  = dI(1) / J_e;
R_e  = (dI(2) - L_e*B_e) / J_e;
Km2  = (dI(3) - R_e*B_e) / Km_e;

fprintf('\n--- MOTOR PARAMETERS ---\n');
fprintf('R=%.6f  L=%.6f  J=%.6f  B=%.6f  Km=%.6f\n', R_e, L_e, J_e, B_e, Km2);

% fine tuning
gTL  = 0.925;
B_f  = 0.295;
TL_s = gTL * load_;

% full simulation with load disturbance
N       = length(time);
ia_s    = zeros(N,1);
wr_s    = zeros(N,1);
ia_s(1) = I0;
wr_s(1) = W0;
for k = 1:N-1
    dia      = -(R_e/L_e)*ia_s(k) - (Km2/L_e)*wr_s(k) + (1/L_e)*volt(k);
    dwr      = (Km_e/J_e)*ia_s(k) - (B_f/J_e)*wr_s(k) - (1/J_e)*TL_s(k);
    ia_s(k+1) = ia_s(k) + ts*dia;
    wr_s(k+1) = wr_s(k) + ts*dwr;
end

fprintf('RMSE speed = %.6f   RMSE current = %.6f\n', ...
    sqrt(mean((wr_s-speed).^2)), sqrt(mean((ia_s-curr).^2)));

% ============================================================
%  FIGURE 1 - Full model vs measured
% ============================================================
c_meas  = [0.85 0.15 0.10];   % rojo medicion
c_model = [0.08 0.45 0.85];   % azul modelo
c_volt  = [0.15 0.15 0.15];
c_load  = [0.75 0.35 0.00];

fig1 = figure('Name','Motor Model vs Measured','Color','w', ...
              'Position',[80 60 920 680]);

ax1 = subplot(3,1,1);
plot(time, speed,  'Color', c_meas,  'LineWidth', 1.6); hold on;
plot(time, wr_s,   'Color', c_model, 'LineWidth', 1.4, 'LineStyle','--');
ylabel('\omega_m  [rad/s]', 'FontSize', 11);
title('Angular Velocity', 'FontSize', 12, 'FontWeight','bold');
legend('Measured', 'Model + disturbance', 'Location','best', 'FontSize', 9);
grid on; box off;
set(ax1, 'GridAlpha', 0.25, 'TickDir','out');

ax2 = subplot(3,1,2);
plot(time, curr,   'Color', c_meas,  'LineWidth', 1.6); hold on;
plot(time, ia_s,   'Color', c_model, 'LineWidth', 1.4, 'LineStyle','--');
ylabel('i_a  [A]', 'FontSize', 11);
title('Armature Current', 'FontSize', 12, 'FontWeight','bold');
legend('Measured', 'Model + disturbance', 'Location','best', 'FontSize', 9);
grid on; box off;
set(ax2, 'GridAlpha', 0.25, 'TickDir','out');

ax3 = subplot(3,1,3);
plot(time, volt,  'Color', c_volt, 'LineWidth', 1.4); hold on;
plot(time, load_, 'Color', c_load, 'LineWidth', 1.4);
ylabel('V_a [V]  /  T_L [Nm]', 'FontSize', 11);
xlabel('Time [s]', 'FontSize', 11);
title('Input Voltage and Load Torque', 'FontSize', 12, 'FontWeight','bold');
legend('V_a', 'T_L', 'Location','best', 'FontSize', 9);
grid on; box off;
set(ax3, 'GridAlpha', 0.25, 'TickDir','out');

% ============================================================
%  FIGURE 2 - Chen identification points
% ============================================================
fig2 = figure('Name','Chen Method - Identification Points','Color','w', ...
              'Position',[100 60 860 520]);

ax4 = subplot(2,1,1);
plot(t_w, spd_w, 'Color',[0.10 0.40 0.75], 'LineWidth', 1.8); hold on;
pts_W = interp1(time, spd_s, resW.tpts_abs, 'linear');
scatter(resW.tpts_abs, pts_W, 80, 'r', 'filled', ...
    'DisplayName','Chen points');
xline(t0, '--', 'Step', 'Color',[0.5 0.5 0.5], 'LabelOrientation','horizontal');
ylabel('\omega_m  [rad/s]', 'FontSize', 11);
title('Speed — Chen Identification Points', 'FontSize', 12, 'FontWeight','bold');
legend('Smoothed speed','Chen points','Location','best','FontSize',9);
grid on; box off;
set(ax4, 'GridAlpha', 0.25, 'TickDir','out');

ax5 = subplot(2,1,2);
plot(t_w, cur_w, 'Color',[0.85 0.45 0.05], 'LineWidth', 1.8); hold on;
pts_I = interp1(time, cur_s, resI.tpts_abs, 'linear');
scatter(resI.tpts_abs, pts_I, 80, 'r', 'filled', ...
    'DisplayName','Chen points');
xline(t0, '--', 'Step', 'Color',[0.5 0.5 0.5], 'LabelOrientation','horizontal');
ylabel('i_a  [A]', 'FontSize', 11);
xlabel('Time [s]', 'FontSize', 11);
title('Current — Chen Identification Points', 'FontSize', 12, 'FontWeight','bold');
legend('Smoothed current','Chen points','Location','best','FontSize',9);
grid on; box off;
set(ax5, 'GridAlpha', 0.25, 'TickDir','out');

% ============================================================
%  FIGURE 3 - Transfer function validation
% ============================================================
fig3 = figure('Name','Transfer Function Validation','Color','w', ...
              'Position',[120 60 860 480]);

ax6 = subplot(2,1,1);
plot(time, speed,  'Color', c_meas,  'LineWidth', 1.6); hold on;
plot(time, spd_tf, 'Color',[0.20 0.65 0.30], 'LineWidth', 1.4, 'LineStyle','-.');
ylabel('\omega_m  [rad/s]', 'FontSize', 11);
title('Speed: Measured vs TF Reconstruction', 'FontSize', 12, 'FontWeight','bold');
legend('Measured','Reconstructed (TF)','Location','best','FontSize',9);
grid on; box off;
set(ax6, 'GridAlpha', 0.25, 'TickDir','out');

ax7 = subplot(2,1,2);
plot(time, curr,   'Color', c_meas,  'LineWidth', 1.6); hold on;
plot(time, cur_tf, 'Color',[0.20 0.65 0.30], 'LineWidth', 1.4, 'LineStyle','-.');
ylabel('i_a  [A]', 'FontSize', 11);
xlabel('Time [s]', 'FontSize', 11);
title('Current: Measured vs TF Reconstruction', 'FontSize', 12, 'FontWeight','bold');
legend('Measured','Reconstructed (TF)','Location','best','FontSize',9);
grid on; box off;
set(ax7, 'GridAlpha', 0.25, 'TickDir','out');

fprintf('\n--- TRANSFER FUNCTIONS ---\n');
disp('Gw(s) ='); disp(Gw)
disp('Gi(s) ='); disp(Gi)
disp('--- ZPK form ---')
disp(zpk(Gw))
disp(zpk(Gi))




% ============================================================
%  LOCAL FUNCTIONS
% ============================================================
function res = findChenPoints(t, y, y0, Ei, K, t0, t_a, t_b, u, ...
                               tmin, tmax, step, tol)
    tc   = tmin:step:tmax;
    best.ok    = false;
    best.nrmse = inf;
    for i = 1:numel(tc)-2
        for j = i+1:numel(tc)-1
            for k = j+1:numel(tc)
                tp = [tc(i), tc(j), tc(k)];
                if tp(1)<t_a || tp(3)>t_b, continue; end
                if abs((tp(2)-tp(1))-(tp(3)-tp(2))) > tol, continue; end
                try
                    tmp = chenFit(t, y, y0, Ei, K, tp, t0, t_a, t_b, u);
                    if isfinite(tmp.nrmse) && tmp.nrmse < best.nrmse
                        best    = tmp;
                        best.ok = true;
                    end
                catch
                end
            end
        end
    end
    if ~best.ok
        error('No valid Chen triplet found.');
    end
    res = best;
end

function res = chenFit(t, y, y0, Ei, K, tp, t0, t_a, t_b, u)
    y1 = interp1(t, y, tp(1),'linear') - y0;
    y2 = interp1(t, y, tp(2),'linear') - y0;
    y3 = interp1(t, y, tp(3),'linear') - y0;

    k1 = y1/(Ei*K) - 1;
    k2 = y2/(Ei*K) - 1;
    k3 = y3/(Ei*K) - 1;

    be = 4*k1^3*k3 - 3*k1^2*k2^2 - 4*k2^3 + k3^2 + 6*k1*k2*k3;
    dn = 2*(k1^2 + k2);
    if be<=0 || abs(dn)<1e-12, error('invalid be'); end

    a1 = (k1*k2 + k3 - sqrt(be)) / dn;
    a2 = (k1*k2 + k3 + sqrt(be)) / dn;
    if ~(isreal(a1)&&isreal(a2)&&a1>0&&a1<1&&a2>0&&a2<1)
        error('invalid alpha');
    end

    beta = (k1 + a2) / (a1 - a2);
    dt   = tp(2) - tp(1);
    T1   = -dt / log(a1);
    T2   = -dt / log(a2);
    T3   = beta*(T1 - T2) + T1;

    s  = tf('s');
    G  = K*(T3*s+1) / ((T1*s+1)*(T2*s+1));

    msk  = (t>=t_a)&(t<=t_b);
    tv   = t(msk);
    yhat = lsim(G, u(msk), tv-tv(1)) + y0;
    yref = y(msk);
    rmse = sqrt(mean((yhat-yref).^2));

    res.tpts_abs = tp;
    res.k1=k1; res.k2=k2; res.k3=k3;
    res.be=be; res.alfa1=a1; res.alfa2=a2;
    res.T1=T1; res.T2=T2; res.T3=T3;
    res.beta=beta; res.G=G; res.rmse=rmse;
    res.nrmse = rmse / max(max(yref)-min(yref), eps);
end