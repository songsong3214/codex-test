function result = No6_RK_fig4()
% No6_RK_fig4
% 用 ode45 按控制方程进行单条联合共振扫频，并保存 MATLAB .fig 文件。
%
% 推荐运行：
%   No4_MCK_w_fig4   % 生成/刷新 No4_coeff_fig4.mat
%   No6_RK_fig4      % 扫频、绘图并保存 Fig4_qualitative_response.fig
%
% 方案A说明：
% - 本文件采用调图用的无量纲模型，不直接使用 No4 中 raw.F1_raw/F2_raw/F3_raw。
% - No4 负责统一模态、保存参数和说明；No6 负责稳健扫频、诊断和出图。
%
% 控制方程：
%   Wddot + theta1*Wdot + omegaL^2*(1+P*cos(Omega*tau))*W
%         + theta2*W^2 + theta3*W^3
%   = kq*cos(0.5*Omega*tau + alpha)
%
% 图中只绘制一条核心曲线：
% combined resonance：P ~= 0 且 kq ~= 0。

clc; close all;

%% 1. 读取 No4 生成的参数；没有 mat 文件时使用同一套默认图4参数
coef = load_or_default_coef();
coef = normalize_coef(coef);
validate_coef(coef);
print_coef_summary(coef);

%% 2. 扫频设置
Omega_arr = build_omega_grid(coef);
setRK = local_rk_settings();

fprintf('\n====== 扫频设置 ======\n');
fprintf('Omega 点数        = %d\n', numel(Omega_arr));
fprintf('Omega 范围        = %.6g 到 %.6g\n', Omega_arr(1), Omega_arr(end));
fprintf('中心频率 2omegaL  = %.6g\n', 2*coef.omegaL);
fprintf('积分周期数        = %d\n', setRK.n_periods);
fprintf('每周期采样点      = %d\n', setRK.points_per_period);
fprintf('稳态截取比例      = %.2f\n', setRK.steady_skip);

%% 3. 核心联合共振曲线：P ~= 0, kq ~= 0
caseCombined = coef;
[Omega_plot, W_combined] = sweep_branch(Omega_arr, [1e-6; 0], caseCombined, setRK, 'forward', 'combined resonance');

%% 4. 诊断和绘图
curveName = 'combinedResonance';
print_curve_summary(Omega_plot, W_combined, curveName);
assert_finite_curve(W_combined, curveName);

fig = figure('Color','w'); hold on; box on;
plot(Omega_plot, W_combined, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Combined resonance');
plot_reference_line(2*coef.omegaL, '2\omega_L');

grid on;
xlabel('\Omega', 'FontSize', 12);
ylabel('W_m', 'FontSize', 12);
title('Qualitative frequency-sweep curve of combined resonance', 'FontSize', 12);
legend('Location', 'best');
xlim([Omega_arr(1), Omega_arr(end)]);
ymax = max(W_combined(:));
if ymax <= 0
    ymax = 1e-4;
end
ylim([0, ymax*1.15 + 1e-5]);

outputFig = 'Fig4_qualitative_response.fig';
save_matlab_figure(fig, outputFig);
fprintf('\n绘图完成，已保存 MATLAB 图文件：%s\n', outputFig);
fprintf('调参建议：峰值位置看 omegaL；幅值大小看 theta3/kq；曲线平滑度看 n_periods。\n');

result = struct();
result.coef = coef;
result.setRK = setRK;
result.Omega = Omega_plot;
result.W_combined = W_combined;
result.figureFile = outputFig;
end

%% ========================================================================


function save_matlab_figure(fig, outputFig)
if exist('savefig', 'file') == 2 || exist('savefig', 'builtin') == 5
    savefig(fig, outputFig);
else
    saveas(fig, outputFig);
end
end

function plot_reference_line(xValue, labelText)
% 兼容较老 MATLAB：新版本用 xline，旧版本退回到普通 plot 竖线。
if exist('xline', 'file') == 2 || exist('xline', 'builtin') == 5
    xline(xValue, 'k:', 'LineWidth', 1.0, 'DisplayName', labelText);
else
    yl = ylim;
    plot([xValue xValue], yl, 'k:', 'LineWidth', 1.0, 'DisplayName', labelText);
    ylim(yl);
end
end

function coef = load_or_default_coef()
if exist('No4_coeff_fig4.mat', 'file') == 2
    S = load('No4_coeff_fig4.mat', 'coef');
    if isfield(S, 'coef')
        coef = S.coef;
        fprintf('已读取 No4_coeff_fig4.mat 中的 coef。\n');
        return;
    end
end
warning('没有找到有效的 No4_coeff_fig4.mat，正在使用 No6 内置默认图4参数。建议先运行 No4_MCK_w_fig4。');
coef = local_default_coef();
end

function coef = normalize_coef(coef)
% 兼容旧 mat 文件：缺字段时补默认值，保证 No6 可以稳定运行。
def = local_default_coef();
fields = fieldnames(def);
for i = 1:numel(fields)
    name = fields{i};
    if ~isfield(coef, name) || isempty(coef.(name))
        coef.(name) = def.(name);
    end
end
coef.F1 = coef.omegaL^2;
coef.Pbar = coef.F1 * coef.P;
coef.F0 = coef.theta1;
coef.R0 = 0;
coef.Fbar = coef.kq;
coef.F2 = coef.theta2;
coef.F3 = coef.theta3;
coef.Omega_center = 2*coef.omegaL;
end

function Omega_arr = build_omega_grid(coef)
if isfield(coef, 'Omega_points') && ~isempty(coef.Omega_points)
    nOmega = coef.Omega_points;
else
    nOmega = 151;
end
Omega_arr = linspace(coef.Omega_min, coef.Omega_max, nOmega);
end

function setRK = local_rk_settings()
setRK.n_periods = 120;
setRK.points_per_period = 60;
setRK.steady_skip = 0.75;
setRK.relTol = 1e-7;
setRK.absTol = 1e-9;
end

function [Omega_plot, Wm_plot] = sweep_branch(Omega_arr, y0, c, setRK, direction, label)
if strcmpi(direction, 'backward')
    Omega_work = fliplr(Omega_arr);
else
    Omega_work = Omega_arr;
end

Wm_work = zeros(size(Omega_work));
y_last = y0(:);
fprintf('开始扫频：%s ...\n', label);

for k = 1:numel(Omega_work)
    Omega = Omega_work(k);
    [Wm_work(k), y_last] = solve_one_frequency(Omega, y_last, c, setRK);
    if mod(k, max(1, floor(numel(Omega_work)/5))) == 0 || k == numel(Omega_work)
        fprintf('  %s: %d/%d, Omega=%.6g, Wm=%.6g\n', label, k, numel(Omega_work), Omega, Wm_work(k));
    end
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
% 因为右端强迫是 cos(Omega*tau/2+alpha)，公共周期取 4*pi/Omega。
T = 4*pi/Omega;
tau_end = setRK.n_periods * T;
N = setRK.n_periods * setRK.points_per_period;
tau_eval = linspace(0, tau_end, N);

opts = odeset('RelTol', setRK.relTol, 'AbsTol', setRK.absTol, 'MaxStep', T/80);
[~, y] = ode45(@(tau,y) reduced_w_ode(tau, y, Omega, c), tau_eval, y0, opts);

if any(~isfinite(y(:)))
    error('Omega=%.6g 时积分出现 NaN/Inf，请减小 kq/P 或增大 theta3。', Omega);
end

idx_start = max(1, floor(setRK.steady_skip * size(y,1)));
W_steady = y(idx_start:end, 1);
Wm = 0.5 * (max(W_steady) - min(W_steady));
y_last = y(end, :).';
end

function dydtau = reduced_w_ode(tau, y, Omega, c)
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
coef.F0 = coef.theta1;
coef.R0 = 0;
coef.omegaL = 0.0192;
coef.F1 = coef.omegaL^2;
coef.P = 1.0;
coef.Pbar = coef.F1 * coef.P;
coef.kq = 0.001;
coef.Fbar = coef.kq;
coef.alpha = -pi/4;
coef.theta2 = 0;
coef.theta3 = 2000;
coef.F2 = coef.theta2;
coef.F3 = coef.theta3;
coef.Omega_center = 2*coef.omegaL;
coef.Omega_span = 0.0024;
coef.Omega_min = coef.Omega_center - coef.Omega_span;
coef.Omega_max = coef.Omega_center + coef.Omega_span;
coef.Omega_points = 151;
coef.model_note = '方案A：定性复现图4的无量纲调参模型。';
end

function validate_coef(c)
need = {'I11','theta1','omegaL','P','kq','alpha','theta2','theta3','Omega_min','Omega_max','Omega_points'};
for i = 1:numel(need)
    if ~isfield(c, need{i})
        error('coef 缺少字段：%s。请重新运行 No4_MCK_w_fig4。', need{i});
    end
    if ~isnumeric(c.(need{i})) || ~isscalar(c.(need{i})) || ~isfinite(c.(need{i}))
        error('coef.%s 必须是有限数值标量。', need{i});
    end
end
if c.I11 <= 0
    error('I11 必须为正。');
end
if c.theta1 < 0
    error('theta1 不能为负。');
end
if c.omegaL <= 0
    error('omegaL 必须为正。');
end
if c.Omega_min >= c.Omega_max
    error('Omega_min 必须小于 Omega_max。');
end
if c.Omega_points < 3 || fix(c.Omega_points) ~= c.Omega_points
    error('Omega_points 必须是不小于 3 的整数。');
end
end

function print_coef_summary(c)
fprintf('\n====== 图4定性复现参数 ======\n');
fprintf('I11          = %.6g\n', c.I11);
fprintf('theta1       = %.6g\n', c.theta1);
fprintf('omegaL       = %.6g\n', c.omegaL);
fprintf('2*omegaL     = %.6g\n', 2*c.omegaL);
fprintf('P            = %.6g\n', c.P);
fprintf('kq           = %.6g\n', c.kq);
fprintf('alpha        = %.6g\n', c.alpha);
fprintf('theta2       = %.6g\n', c.theta2);
fprintf('theta3       = %.6g\n', c.theta3);
fprintf('Omega_min/max= %.6g / %.6g\n', c.Omega_min, c.Omega_max);
end

function print_curve_summary(Omega_arr, Wm, curveName)
fprintf('\n====== 曲线诊断 ======\n');
[mx, idx] = max(Wm);
fprintf('%-24s max Wm = %.6g at Omega = %.6g\n', curveName, mx, Omega_arr(idx));
end

function assert_finite_curve(Wm, curveName)
if any(~isfinite(Wm(:)))
    error('%s 曲线包含 NaN/Inf。', curveName);
end
end
