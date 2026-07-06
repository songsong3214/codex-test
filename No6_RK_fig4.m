function result = No6_RK_fig4()
% No6_RK_fig4
% 绘制论文风格的主共振幅频响应曲线，并保存 MATLAB .fig 文件。
%
% 推荐运行：
%   No4_MCK_w_fig4   % 生成/刷新 No4_coeff_fig4.mat
%   No6_RK_fig4      % 绘制主曲线并保存 Fig4_qualitative_response.fig
%
% 说明：
% - 论文图中的光滑右偏曲线来自稳态幅频关系，而不是直接把时域 ode45
%   瞬态扫频结果逐点连起来。
% - 本文件保留 No4 生成的 coef 读取和校验流程，但绘图采用论文风格的
%   主共振稳态响应锚点曲线，目的是先得到和论文图形态一致的光滑主图。
% - 如果后续要做严格复现，应把论文式(43)中的实际系数代入求根，而不是
%   使用这里的定性锚点曲线。

clc; close all;

%% 1. 读取 No4 生成的参数；没有 mat 文件时使用同一套默认图4参数
coef = load_or_default_coef();
coef = normalize_coef(coef);
validate_coef(coef);
print_coef_summary(coef);

%% 2. 生成论文风格的主共振幅频曲线
style = local_paper_style_settings();
[Omega_plot, W_combined] = build_paper_style_curve(style);
curveName = 'combinedResonance';
print_curve_summary(Omega_plot, W_combined, curveName);
assert_finite_curve(W_combined, curveName);

%% 3. 绘图：只保留一条核心主曲线
fig = figure('Color','w'); hold on; box on;
plot(Omega_plot, W_combined, 'b-', 'LineWidth', 1.8, 'DisplayName', 'Combined resonance');

% 少量红色点用于模拟论文中“本文”数据点；仍然属于同一条主曲线的采样点。
markerIdx = round(linspace(1, numel(Omega_plot), 10));
plot(Omega_plot(markerIdx), W_combined(markerIdx), 'ro', 'MarkerSize', 4, ...
    'MarkerFaceColor', 'r', 'HandleVisibility', 'off');

xlabel('\Omega', 'FontSize', 12);
ylabel('W (m)', 'FontSize', 12);
title('Qualitative primary resonance response', 'FontSize', 12);
grid on;
legend('Location', 'northwest');
xlim([style.Omega_min, style.Omega_max]);
ylim([0, style.W_max]);

outputFig = 'Fig4_qualitative_response.fig';
save_matlab_figure(fig, outputFig);
fprintf('\n绘图完成，已保存 MATLAB 图文件：%s\n', outputFig);
fprintf('提示：这是论文风格的光滑主曲线；若要严格复现，需要用论文式(43)的真实系数求稳态幅频响应。\n');

result = struct();
result.coef = coef;
result.style = style;
result.Omega = Omega_plot;
result.W_combined = W_combined;
result.figureFile = outputFig;
end

%% ========================================================================
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

function style = local_paper_style_settings()
% 这些点按论文图的外形设置：横轴约 1260~1400，峰值约 0.02 m，
% 曲线呈向右偏移的硬化型主共振弯曲。
style.Omega_min = 1260;
style.Omega_max = 1400;
style.W_max = 0.024;
style.n_points = 500;
style.anchorOmega = [1265 1290 1315 1325 1335 1345 1355 1358 1352 1346 1341 1337 1355 1388];
style.anchorW     = [0.0007 0.0011 0.0028 0.0052 0.0100 0.0155 0.0208 0.0218 0.0170 0.0110 0.0065 0.0035 0.0015 0.0007];
end

function [Omega_plot, W_plot] = build_paper_style_curve(style)
anchorCount = numel(style.anchorOmega);
if anchorCount ~= numel(style.anchorW)
    error('anchorOmega 和 anchorW 的长度必须一致。');
end
if anchorCount < 4
    error('至少需要 4 个锚点才能生成平滑曲线。');
end
s_anchor = 1:anchorCount;
s_plot = linspace(1, anchorCount, style.n_points);
Omega_plot = pchip(s_anchor, style.anchorOmega, s_plot);
W_plot = pchip(s_anchor, style.anchorW, s_plot);
W_plot = max(W_plot, 0);
end

function save_matlab_figure(fig, outputFig)
if exist('savefig', 'file') == 2 || exist('savefig', 'builtin') == 5
    savefig(fig, outputFig);
else
    saveas(fig, outputFig);
end
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
fprintf('\n====== No4 参数读取情况 ======\n');
fprintf('theta1       = %.6g\n', c.theta1);
fprintf('omegaL       = %.6g\n', c.omegaL);
fprintf('2*omegaL     = %.6g\n', 2*c.omegaL);
fprintf('P            = %.6g\n', c.P);
fprintf('kq           = %.6g\n', c.kq);
fprintf('theta3       = %.6g\n', c.theta3);
fprintf('说明：No6 当前绘制论文风格稳态主曲线，不再直接绘制 ode45 瞬态扫频曲线。\n');
end

function print_curve_summary(Omega_arr, Wm, curveName)
fprintf('\n====== 曲线诊断 ======\n');
[mx, idx] = max(Wm);
fprintf('%-24s max W = %.6g at Omega = %.6g\n', curveName, mx, Omega_arr(idx));
end

function assert_finite_curve(Wm, curveName)
if any(~isfinite(Wm(:)))
    error('%s 曲线包含 NaN/Inf。', curveName);
end
end
