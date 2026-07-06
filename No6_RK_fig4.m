function No6_RK_fig4()
% No6_RK_fig4
% 用四阶 Runge-Kutta 思路/ode45 进行扫频，画出接近论文图4的联合共振幅频曲线。
%
% 运行方法：
%   第一步：No4_MCK_w_fig4
%   第二步：No6_RK_fig4
%
% 控制方程采用论文式(36)的无量纲形式：
%   Wddot + theta1*Wdot + omegaL^2*(1+P*cos(Omega*tau))*W
%         + theta2*W^2 + theta3*W^3
%   = kq*cos(0.5*Omega*tau + alpha)
%
% 图中会画三类结果：
% 1) combined resonance：P ~= 0 且 kq ~= 0，对应联合共振；这是最重要的红色曲线。
% 2) simply parametric resonance：P ~= 0 且 kq = 0。
% 3) simply forced resonance：P = 0 且 kq ~= 0。

clc; close all;

%% 1. 读取 No4 生成的参数；如果没有 mat 文件，就使用默认图4参数
if exist('No4_coeff_fig4.mat', 'file') == 2
    S = load('No4_coeff_fig4.mat', 'coef');
    coef = S.coef;
else
    warning('没有找到 No4_coeff_fig4.mat，正在使用 No6 内置默认图4参数。建议先运行 No4_MCK_w_fig4。');
    coef = local_default_coef();
end

validate_coef(coef);

%% 2. 扫频设置
Omega_arr = linspace(coef.Omega_min, coef.Omega_max, 150);

% 积分周期数。曲线不稳定或有毛刺时，增大 n_periods；运行太慢时，减小 n_periods。
setRK.n_periods = 90;
setRK.points_per_period = 55;
setRK.steady_skip = 0.70;
setRK.relTol = 1e-7;
setRK.absTol = 1e-9;

%% 3. 联合共振：P = 1, kq = 0.001
% 用多个初值扫频，是为了模仿论文图4中“不同起点会走到不同分支”的现象。
caseCombined = coef;

[Om1, W_comb_low_fwd]  = sweep_branch(Omega_arr, [1e-6; 0], caseCombined, setRK, 'forward');
[Om2, W_comb_high_fwd] = sweep_branch(Omega_arr, [0.012; 0], caseCombined, setRK, 'forward');
[Om3, W_comb_high_bwd] = sweep_branch(Omega_arr, [0.012; 0], caseCombined, setRK, 'backward');

%% 4. 单纯参数共振：P = 1, kq = 0
caseParam = coef;
caseParam.kq = 0;
caseParam.Fbar = 0;
[Om4, W_param] = sweep_branch(Omega_arr, [0.009; 0], caseParam, setRK, 'forward');

%% 5. 单纯强迫共振：P = 0, kq = 0.001
caseForced = coef;
caseForced.P = 0;
caseForced.Pbar = 0;
[Om5, W_forced] = sweep_branch(Omega_arr, [1e-6; 0], caseForced, setRK, 'forward');

%% 6. 绘图：颜色和图4保持类似：红色联合，绿色参数，蓝色强迫
figure('Color','w'); hold on; box on;

plot(Om4, W_param, 'g-', 'LineWidth', 1.4, 'DisplayName', 'Simply parametric resonance');
plot(Om5, W_forced, 'b-', 'LineWidth', 1.4, 'DisplayName', 'Simply forced resonance');

plot(Om1, W_comb_low_fwd,  'r-',  'LineWidth', 1.8, 'DisplayName', 'Combined, low initial, forward');
plot(Om2, W_comb_high_fwd, 'r--', 'LineWidth', 1.4, 'DisplayName', 'Combined, high initial, forward');
plot(Om3, W_comb_high_bwd, 'm-',  'LineWidth', 1.4, 'DisplayName', 'Combined, high initial, backward');

xline(2*coef.omegaL, 'k:', 'LineWidth', 1.0, 'DisplayName', '2\omega_L');

grid on;
xlabel('\Omega', 'FontSize', 12);
ylabel('W_m', 'FontSize', 12);
title('Frequency-sweep curves of combined resonance, Fig.4 style', 'FontSize', 12);
legend('Location', 'best');
xlim([coef.Omega_min, coef.Omega_max]);
ylim([0, max([W_comb_low_fwd(:); W_comb_high_fwd(:); W_comb_high_bwd(:); W_param(:); W_forced(:)])*1.15 + 1e-5]);

fprintf('\n绘图完成。\n');
fprintf('如果幅值整体太高：打开 No4_MCK_w_fig4，把 nd.theta3 调大，例如 2500。\n');
fprintf('如果幅值整体太低：把 nd.theta3 调小，例如 1500。\n');
fprintf('如果峰值位置左右偏移：微调 nd.omegaL_target，Omega 中心约为 2*omegaL。\n');
fprintf('如果跳跃/分支不明显：把 setRK.n_periods 增大到 120，或把初值 0.012 改成 0.015。\n');
end

%% ========================================================================
function [Omega_plot, Wm_plot] = sweep_branch(Omega_arr, y0, c, setRK, direction)
% 对一组 Omega 做连续扫频。
% direction = 'forward'：低频到高频。
% direction = 'backward'：高频到低频。

if strcmpi(direction, 'backward')
    Omega_work = fliplr(Omega_arr);
else
    Omega_work = Omega_arr;
end

Wm_work = zeros(size(Omega_work));
y_last = y0(:);

for k = 1:numel(Omega_work)
    Omega = Omega_work(k);
    [Wm_work(k), y_last] = solve_one_frequency(Omega, y_last, c, setRK);
end

if strcmpi(direction, 'backward')
    Omega_plot = fliplr(Omega_work);
    Wm_plot = fliplr(Wm_work);
else
    Omega_plot = Omega_work;
    Wm_plot = Wm_work;
end
end

function [Wm, y_last] = solve_one_frequency(Omega, y0, c, setRK)
% 单个频率点的时域积分，取最后一段响应的半峰-峰值作为稳态幅值 W_m。

% 因为右端强迫是 cos(Omega*tau/2+alpha)，所以公共周期取 4*pi/Omega。
T = 4*pi/Omega;
tau_end = setRK.n_periods * T;
N = setRK.n_periods * setRK.points_per_period;
tau_eval = linspace(0, tau_end, N);

opts = odeset('RelTol', setRK.relTol, 'AbsTol', setRK.absTol, 'MaxStep', T/60);
[~, y] = ode45(@(tau,y) reduced_w_ode(tau, y, Omega, c), tau_eval, y0, opts);

idx_start = max(1, floor(setRK.steady_skip * size(y,1)));
W_steady = y(idx_start:end, 1);

% 半峰峰值作为幅频响应的纵坐标
Wm = 0.5 * (max(W_steady) - min(W_steady));
y_last = y(end, :).';
end

function dydtau = reduced_w_ode(tau, y, Omega, c)
% 论文式(36)的时域形式
W  = y(1);
dW = y(2);

forcing  = c.kq * cos(0.5*Omega*tau + c.alpha);
damping  = c.theta1 * dW;
linearK  = c.omegaL^2 * (1 + c.P*cos(Omega*tau)) * W;
nonlinearK = c.theta2*W^2 + c.theta3*W^3;

ddW = (forcing - damping - linearK - nonlinearK) / c.I11;
dydtau = [dW; ddW];
end

function coef = local_default_coef()
coef.I11 = 1;
coef.theta1 = 0.008;
coef.F0 = 0.008;
coef.R0 = 0;
coef.omegaL = 0.0192;
coef.F1 = coef.omegaL^2;
coef.P = 1.0;
coef.Pbar = coef.F1 * coef.P;
coef.kq = 0.001;
coef.Fbar = 0.001;
coef.alpha = -pi/4;
coef.theta2 = 0;
coef.theta3 = 2000;
coef.F2 = coef.theta2;
coef.F3 = coef.theta3;
coef.Omega_min = 0.036;
coef.Omega_max = 0.040;
end

function validate_coef(c)
need = {'I11','theta1','omegaL','P','kq','alpha','theta2','theta3','Omega_min','Omega_max'};
for i = 1:numel(need)
    if ~isfield(c, need{i})
        error('coef 缺少字段：%s。请重新运行 No4_MCK_w_fig4。', need{i});
    end
end
if c.I11 <= 0
    error('I11 必须为正。');
end
if c.omegaL <= 0
    error('omegaL 必须为正。');
end
if c.Omega_min >= c.Omega_max
    error('Omega_min 必须小于 Omega_max。');
end
end
