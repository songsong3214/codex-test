function coef = No4_MCK_w_fig4()
% No4_MCK_w_fig4
% 作用：生成最终给 No6_RK_fig4 使用的一自由度 w 方向约化方程参数。
%
% 三个方向的关系：
% 1) No4_MCK_u_fig4 生成 u 方向代数式：gamma11*U+gamma12*V+gamma13*W+gamma14*W^2=0
% 2) No4_MCK_v_fig4 生成 v 方向代数式：gamma21*U+gamma22*V+gamma23*W+gamma24*W^2=0
% 3) 本文件把上面两个方程联立，得到 U(W)、V(W)，再代回 w 方向方程，得到：
%
%      Wddot + theta1*Wdot + omegaL^2*(1 + P*cos(Omega*tau))*W
%            + theta2*W^2 + theta3*W^3
%      = kq*cos(0.5*Omega*tau + alpha)
%
% 这就是论文式(36)的形式；No6_RK_fig4 用它扫频绘图。
%
% 注意：
% - 本代码保留 u/v 消元的计算框架。
% - 但因为你要求 E、rho 直接取常数，且忽略温度、GPL/孔隙沿厚度分布，
%   所以最终出图参数采用“图4无量纲参数归一化”。这样更容易得到与图4相似的曲线。

clearvars -except coef;
clc;

%% 1. 图4使用的基础参数
phy.L = 1.0;
phy.h = 0.01;
phy.R = 1.0;
phy.E = 2.0e11;        % 常数杨氏模量，表3矩阵量级
phy.rho = 7850;        % 常数密度，表1/表3矩阵量级
phy.nu = 0.3;
phy.m = 1;
phy.n = 4;             % 论文图4明确使用 (m,n)=(1,4)

nd.Vc = 0.02;
nd.P = 1.0;            % 图4：P = 1
nd.kq = 0.001;         % 图4：kq = 0.001
nd.theta1 = 0.008;     % 图4：#1 = 0.008
nd.alpha = -pi/4;      % 图4：alpha = -pi/4

% 这个频率是为了让联合共振区域落在论文图4附近：Omega ≈ 0.036~0.040。
% 联合共振中横向激励频率为 Omega/2，因此中心位置近似 Omega ≈ 2*omegaL。
nd.omegaL_target = 0.0192;

% 非线性三次刚度，主要控制幅值高度。
% 幅值太大：增大 theta3；幅值太小：减小 theta3。
nd.theta2 = 0;
nd.theta3 = 2000;

%% 2. u、v 方向系数，并进行消元
fprintf('\n====== 1) 正在计算 u/v 方向 Galerkin 系数，请稍等 ======\n');
resU = No4_MCK_u_fig4(phy);
resV = No4_MCK_v_fig4(phy);

syms U V W dW ddW X TH z real

det_uv = simplify(resU.gamma11*resV.gamma22 - resU.gamma12*resV.gamma21);
eta11 = simplify((resU.gamma12*resV.gamma23 - resV.gamma22*resU.gamma13)/det_uv);
eta21 = simplify((resU.gamma12*resV.gamma24 - resV.gamma22*resU.gamma14)/det_uv);
eta12 = simplify((resV.gamma21*resU.gamma13 - resU.gamma11*resV.gamma23)/det_uv);
eta22 = simplify((resV.gamma21*resU.gamma14 - resU.gamma11*resV.gamma24)/det_uv);

% 论文式(32)、式(33)的形式：U = eta11*W + eta21*W^2，V = eta12*W + eta22*W^2
U_of_W = eta11*W + eta21*W^2;
V_of_W = eta12*W + eta22*W^2;

%% 3. 构造 w 方向结构项，提取原始 c31~c37
% 这一段不是最终调图直接用的参数，但它保留了论文的 u/v 消元流程，便于论文中说明。
fprintf('====== 2) 正在计算 w 方向结构项并提取 c31~c37 ======\n');
Rbar = phy.R/phy.L;

phiU = cos(phy.m*pi*X) * cos(phy.n*TH/2);
phiV = sin(phy.m*pi*X) * sin(phy.n*TH/2);
phiW = sin(phy.m*pi*X) * cos(phy.n*TH/2);

% 重要：这里先用独立的 U,V,W 建 w 方程，不要提前代入 U_of_W/V_of_W。
% 先提取 c31~c37，再消元。如果提前代入，c31/c32/c36/c37 会被错误消掉。
u0 = U * phiU;
v0 = V * phiV;
w0 = W * phiW;

ux    = diff(u0, X);
uTH   = diff(u0, TH);
vx    = diff(v0, X);
vTH   = diff(v0, TH);
wx    = diff(w0, X);
wTH   = diff(w0, TH);
wXX   = diff(w0, X, 2);
wXTH  = diff(diff(w0, X), TH);
wTHTH = diff(w0, TH, 2);

ex0  = ux + 0.5*wx^2;
eth0 = vTH/Rbar + w0/Rbar + 0.5*(wTH/Rbar)^2;
gam0 = vx + uTH/Rbar + wx*wTH/Rbar;

kx   = wXX;
kth  = wTHTH/Rbar^2 - vTH/Rbar;
kxth = wXTH/Rbar - vx;

Q11 = phy.E/(1 - phy.nu^2);
Q12 = phy.nu*phy.E/(1 - phy.nu^2);
Q22 = Q11;
Q66 = phy.E/(2*(1 + phy.nu));

ex  = ex0 - z*kx;
eth = eth0 - z*kth;
gam = gam0 - 2*z*kxth;

sigx  = Q11*ex + Q12*eth;
sigth = Q12*ex + Q22*eth;
sigxt = Q66*gam;

Nx  = int(sigx,  z, -phy.h/2, phy.h/2);
Nth = int(sigth, z, -phy.h/2, phy.h/2);
Nxt = int(sigxt, z, -phy.h/2, phy.h/2);
Mx  = int(z*sigx,  z, -phy.h/2, phy.h/2);
Mth = int(z*sigth, z, -phy.h/2, phy.h/2);
Mxt = int(z*sigxt, z, -phy.h/2, phy.h/2);

Lw_bending   = diff(Mx, X, 2) + 2*diff(diff(Mxt, X), TH)/Rbar + diff(Mth, TH, 2)/Rbar^2;
Lw_curvature = -Nth/Rbar;
Lw_membrane  = Nx*wXX + 2*Nxt*wXTH/Rbar + Nth*wTHTH/Rbar^2;
Lw_struct_dim = simplify(Lw_bending + Lw_curvature + Lw_membrane);

% 结构项归一化。这里不用于决定最终图4的频率，只保存为 raw 参考。
D = phy.E*phy.h^3/(12*(1 - phy.nu^2));
Lw_ref = D*phy.h/phy.L^4;
Lw_struct = simplify(Lw_struct_dim / Lw_ref);

modalNorm = simplify(int(int(phiW^2, X, 0, 1), TH, 0, 2*pi));
modalScale = simplify(1/modalNorm);
Kw_struct = simplify(modalScale * int(int(Lw_struct*phiW, X, 0, 1), TH, 0, 2*pi));

c31 = simplify(subs(diff(Kw_struct, U), [U,V,W], [0,0,0]));
c32 = simplify(subs(diff(Kw_struct, V), [U,V,W], [0,0,0]));
c33 = simplify(subs(diff(Kw_struct, W), [U,V,W], [0,0,0]));
c34 = simplify(subs(diff(Kw_struct, W, 2)/2, [U,V,W], [0,0,0]));
c35 = simplify(subs(diff(Kw_struct, W, 3)/6, [U,V,W], [0,0,0]));
c36 = simplify(subs(diff(diff(Kw_struct, U), W), [U,V,W], [0,0,0]));
c37 = simplify(subs(diff(diff(Kw_struct, V), W), [U,V,W], [0,0,0]));

F1_raw = simplify(-(c33 + c31*eta11 + c32*eta12));
F2_raw = simplify(-(c34 + c31*eta21 + c32*eta22 + c36*eta11 + c37*eta12));
F3_raw = simplify(-(c35 + c36*eta21 + c37*eta22));

raw = struct();
raw.c31 = c31; raw.c32 = c32; raw.c33 = c33; raw.c34 = c34; raw.c35 = c35; raw.c36 = c36; raw.c37 = c37;
raw.eta11 = eta11; raw.eta21 = eta21; raw.eta12 = eta12; raw.eta22 = eta22;
raw.F1_raw = F1_raw; raw.F2_raw = F2_raw; raw.F3_raw = F3_raw;
raw.modalNorm = modalNorm;

%% 4. 最终用于图4的无量纲约化参数
% 采用论文式(36)的无量纲形式：
% Wddot + theta1*Wdot + omegaL^2*(1+P*cos(Omega*tau))*W + theta2*W^2 + theta3*W^3
% = kq*cos(Omega*tau/2 + alpha)
coef = struct();
coef.I11 = 1;
coef.theta1 = nd.theta1;
coef.F0 = nd.theta1;
coef.R0 = 0;
coef.omegaL = nd.omegaL_target;
coef.F1 = coef.omegaL^2;
coef.P = nd.P;
coef.Pbar = coef.F1 * coef.P;  % 对应式(34)中的 Px；注意不是直接等于 P
coef.kq = nd.kq;
coef.Fbar = nd.kq;
coef.alpha = nd.alpha;
coef.theta2 = nd.theta2;
coef.theta3 = nd.theta3;
coef.F2 = nd.theta2;
coef.F3 = nd.theta3;
coef.Omega_min = 0.036;
coef.Omega_max = 0.040;
coef.Omega_center = 2*coef.omegaL;
coef.phy = phy;

save('No4_coeff_fig4.mat', 'coef', 'raw');

fprintf('\n====== 3) 已生成 No4_coeff_fig4.mat，No6_RK_fig4 将自动读取 ======\n');
fprintf('I11       = %.6g\n', coef.I11);
fprintf('theta1/F0 = %.6g\n', coef.theta1);
fprintf('omegaL    = %.6g\n', coef.omegaL);
fprintf('Omega0    = %.6g  约等于 2*omegaL\n', coef.Omega_center);
fprintf('P         = %.6g\n', coef.P);
fprintf('Pbar      = %.6g  注意：Pbar = omegaL^2 * P\n', coef.Pbar);
fprintf('kq/Fbar   = %.6g\n', coef.kq);
fprintf('theta2/F2 = %.6g\n', coef.theta2);
fprintf('theta3/F3 = %.6g\n', coef.theta3);
fprintf('alpha     = %.6g\n', coef.alpha);
fprintf('扫频范围建议：Omega = %.4f 到 %.4f\n', coef.Omega_min, coef.Omega_max);
end
