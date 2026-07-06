function out = No4_MCK_v_fig4(phy)
% No4_MCK_v_fig4
% v 方向 Galerkin 准静态方程：
%   gamma21*U + gamma22*V + gamma23*W + gamma24*W^2 = 0
%
% 说明：
% 1) 本文件对应论文中式(26)的 V 模态，以及式(29)的 v 方向代数方程。
% 2) 材料按常数处理；不考虑温度项。
% 3) 本文件只负责生成 v 方程的 4 个 gamma 系数，供 No4_MCK_w_fig4 调用。

if nargin < 1 || isempty(phy)
    phy = local_default_phy();
end

Rbar = phy.R / phy.L;

syms X TH U V W z real

phiU = cos(phy.m*pi*X) * cos(phy.n*TH/2);
phiV = sin(phy.m*pi*X) * sin(phy.n*TH/2);
phiW = sin(phy.m*pi*X) * cos(phy.n*TH/2);

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
wTHTH = diff(w0, TH, 2);
wxTH  = diff(diff(w0, X), TH);

ex0  = ux + 0.5*wx^2;
eth0 = vTH/Rbar + w0/Rbar + 0.5*(wTH/Rbar)^2;
gam0 = vx + uTH/Rbar + wx*wTH/Rbar;

kx   = wXX;
kth  = wTHTH/Rbar^2 - vTH/Rbar;
kxth = wxTH/Rbar - vx;

Q11 = phy.E / (1 - phy.nu^2);
Q12 = phy.nu * phy.E / (1 - phy.nu^2);
Q22 = Q11;
Q66 = phy.E / (2*(1 + phy.nu));

ex  = ex0 - z*kx;
eth = eth0 - z*kth;
gam = gam0 - 2*z*kxth;

sigth = Q12*ex + Q22*eth;
sigxt = Q66*gam;

Nth = int(sigth, z, -phy.h/2, phy.h/2);
Nxt = int(sigxt, z, -phy.h/2, phy.h/2);
Mth = int(z*sigth, z, -phy.h/2, phy.h/2);
Mxt = int(z*sigxt, z, -phy.h/2, phy.h/2);

% v 方向平衡方程，对应论文式(19)左端
Rv_membrane = diff(Nxt, X) + diff(Nth, TH)/Rbar;
Rv_bending  = (2*diff(Mxt, X) + diff(Mth, TH)/Rbar) / Rbar;
Rv = simplify(Rv_membrane + Rv_bending);

GalerkinV = int(int(Rv * phiV, X, 0, 1), TH, 0, 2*pi);
GalerkinV = simplify(GalerkinV);

% 提取式(29)的系数
out.gamma21 = simplify(subs(diff(GalerkinV, U), [U,V,W], [0,0,0]));
out.gamma22 = simplify(subs(diff(GalerkinV, V), [U,V,W], [0,0,0]));
out.gamma23 = simplify(subs(diff(GalerkinV, W), [U,V,W], [0,0,0]));
out.gamma24 = simplify(subs(diff(GalerkinV, W, 2)/2, [U,V,W], [0,0,0]));
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
phy.n = 4;
end
