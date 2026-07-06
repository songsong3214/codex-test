function out = No4_MCK_u_fig4(phy)
% No4_MCK_u_fig4
% u 方向 Galerkin 准静态方程：
%   gamma11*U + gamma12*V + gamma13*W + gamma14*W^2 = 0
%
% 说明：
% 1) 本文件对应论文中式(25)的 U 模态，以及式(28)的 u 方向代数方程。
% 2) 材料按你的要求简化为常数 E、rho、nu；不考虑温度项。
% 3) 本文件只负责生成 u 方程的 4 个 gamma 系数，供 No4_MCK_w_fig4 调用。

if nargin < 1 || isempty(phy)
    phy = local_default_phy();
end

Rbar = phy.R / phy.L;      % 无量纲半径 R/L

syms X TH U V W z real

% 论文式(25)~(27)：注意环向模态是 cos(n*TH/2)、sin(n*TH/2)
phiU = cos(phy.m*pi*X) * cos(phy.n*TH/2);
phiV = sin(phy.m*pi*X) * sin(phy.n*TH/2);
phiW = sin(phy.m*pi*X) * cos(phy.n*TH/2);

u0 = U * phiU;
v0 = V * phiV;
w0 = W * phiW;

% Donnell 非线性应变，形式对应论文式(9)、式(10)
ux    = diff(u0, X);
uTH   = diff(u0, TH);
vx    = diff(v0, X);
vTH   = diff(v0, TH);
wx    = diff(w0, X);
wTH   = diff(w0, TH);
wXX   = diff(w0, X, 2);
wTHTH = diff(w0, TH, 2);
wxTH  = diff(diff(w0, X), TH);

ex0  = ux + 0.5*wx^2;
eth0 = vTH/Rbar + w0/Rbar + 0.5*(wTH/Rbar)^2;
gam0 = vx + uTH/Rbar + wx*wTH/Rbar;

kx   = wXX;
kth  = wTHTH/Rbar^2 - vTH/Rbar;
kxth = wxTH/Rbar - vx;

% 各向同性常材料刚度；忽略论文式(11)中的温度项
Q11 = phy.E / (1 - phy.nu^2);
Q12 = phy.nu * phy.E / (1 - phy.nu^2);
Q66 = phy.E / (2*(1 + phy.nu));

ex  = ex0 - z*kx;
eth = eth0 - z*kth;
gam = gam0 - 2*z*kxth;

sigx  = Q11*ex + Q12*eth;
sigxt = Q66*gam;

Nx  = int(sigx,  z, -phy.h/2, phy.h/2);
Nxt = int(sigxt, z, -phy.h/2, phy.h/2);

% u 方向平衡方程，对应论文式(18)左端；u/v 在这里按准静态处理
Ru = diff(Nx, X) + diff(Nxt, TH)/Rbar;

% Galerkin 投影：X 从 0 到 1，TH 从 0 到 2*pi
GalerkinU = int(int(Ru * phiU, X, 0, 1), TH, 0, 2*pi);
GalerkinU = simplify(GalerkinU);

% 提取式(28)的系数
out.gamma11 = simplify(subs(diff(GalerkinU, U), [U,V,W], [0,0,0]));
out.gamma12 = simplify(subs(diff(GalerkinU, V), [U,V,W], [0,0,0]));
out.gamma13 = simplify(subs(diff(GalerkinU, W), [U,V,W], [0,0,0]));
out.gamma14 = simplify(subs(diff(GalerkinU, W, 2)/2, [U,V,W], [0,0,0]));
out.phy = phy;
end

function phy = local_default_phy()
phy.L = 1.0;
phy.h = 0.01;
phy.R = 1.0;
phy.E = 2.0e11;
phy.rho = 7850;
phy.nu = 0.3;
phy.m = 1;
phy.n = 4;       % 图4使用 (m,n)=(1,4)
end
