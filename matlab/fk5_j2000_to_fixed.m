function [C, info] = fk5_j2000_to_fixed(utc, dat, dut1, pm, dn)
% IAU 1976 precession, full 106-term IAU 1980 nutation, IERS 1996 eqeq.
% nut80.dat source: CelesTrak/fundamentals-of-astrodynamics/datalib.
% No Aerospace Toolbox or JPL ephemeris data package is required.
%
% 输入:
%   utc  : UTC 历元 [Y M D h m s]
%   dat  : TAI-UTC 跳秒 (s)
%   dut1 : UT1-UTC (s)
%   pm   : 极移 [xp yp] (rad)
%   dn   : 额外章动修正 [dPsi dEpsilon] (rad)
% 输出:
%   C    : J2000 -> 地固系(ITRF) 的旋转矩阵
%   info : 中间量(时间参数、分矩阵等)

% ---------- 1. 时间参数 ----------
jd0 = juliandate(utc(1), utc(2), utc(3));               % 当日 0h 的儒略日
s   = utc(4) * 3600 + utc(5) * 60 + utc(6);             % 当日已过的秒数

% TT 世纪数：UTC + dat + 32.184 = TT，用于岁差与章动
T  = ((jd0 - 2451545) + (s + dat + 32.184) / 86400) / 36525;

% UT1 世纪数：UTC + dut1 = UT1，用于地球自转(GMST)
Tu = ((jd0 - 2451545) + (s + dut1) / 86400) / 36525;

as2r = pi / 648000;   % 角秒 -> 弧度

% ---------- 2. 旋转矩阵辅助函数 ----------
Rz = @(x) [cos(x) -sin(x) 0; sin(x) cos(x) 0; 0 0 1];
Ry = @(x) [cos(x) 0 sin(x); 0 1 0; -sin(x) 0 cos(x)];
Rx = @(x) [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];

% ---------- 3. 岁差矩阵 P (IAU 1976) ----------
zeta  = (2306.2181 * T + 0.30188 * T^2 + 0.017998 * T^3) * as2r;
theta = (2004.3109 * T - 0.42665 * T^2 - 0.041833 * T^3) * as2r;
z     = (2306.2181 * T + 1.09468 * T^2 + 0.018203 * T^3) * as2r;
P = Rz(z) * Ry(-theta) * Rz(zeta);   % 3-2-3 转序: J2000 平赤道 -> MOD

% ---------- 4. 章动矩阵 N (IAU 1980, 106 项) ----------
% 5 个 Delaunay 基本幅角 (单位: 角秒 -> 弧度)
args = [134.96298139 * 3600 + 1717915922.633 * T + 31.310 * T^2 + 0.064 * T^3; ...  % l  月球平近点角
        357.52772333 * 3600 + 129596581.224 * T -  0.577 * T^2 - 0.012 * T^3; ...  % l' 太阳平近点角
         93.27191028 * 3600 + 1739527263.137 * T - 13.257 * T^2 + 0.011 * T^3; ...  % F  月球纬度幅角
        297.85036306 * 3600 + 1602961601.328 * T -  6.891 * T^2 + 0.019 * T^3; ...  % D  日月平距
        125.04452222 * 3600 - 6962890.539 * T +  7.455 * T^2 + 0.008 * T^3] * as2r; % Om 月球升交点黄经

% 读取 106 项章动表(nut80.dat)，用 persistent 只读一次盘
persistent tab
if isempty(tab)
    tab = load(fullfile(fileparts(mfilename('fullpath')), 'nut80.dat'));
    if size(tab, 1) ~= 106 || size(tab, 2) < 9
        error('Alignment:NutationTable', 'Expected full 106-term IAU 1980 table.');
    end
end

% 每项相位 = 5 个幅角的线性组合
phase = tab(:, 1:5) * args;

% 黄经章动 dPsi 与交角章动 dEpsilon (表内单位 1e-4 角秒)
dpsi = sum((tab(:, 6) + tab(:, 7) * T) .* sin(phase)) * 1e-4 * as2r + dn(1);
deps = sum((tab(:, 8) + tab(:, 9) * T) .* cos(phase)) * 1e-4 * as2r + dn(2);

% 平黄赤交角 eps (IAU 1976)
eps = (84381.448 - 46.8150 * T - 0.00059 * T^2 + 0.001813 * T^3) * as2r;

% 章动矩阵: 3-1-3 转序, MOD -> TOD
N = Rx(eps + deps) * Rz(dpsi) * Rx(-eps);

% ---------- 5. 地球自转: GMST / GAST ----------
% 格林尼治平恒星时 GMST (IAU 1982)，秒 -> 弧度
gmst = mod((67310.54841 + (876600 * 3600 + 8640184.812866) * Tu ...
           + 0.093104 * Tu^2 - 6.2e-6 * Tu^3) * pi / 43200, 2 * pi);

% 春分点方程(equation of the equinoxes)
eq = dpsi * cos(eps);

% IERS 1996: 1997-01-01 之后追加两个补充分项
if jd0 + (s + dut1) / 86400 > 2450449.5
    eq = eq + (0.00264 * sin(args(5)) + 0.000063 * sin(2 * args(5))) * as2r;
end

% 格林尼治真恒星时 GAST
gast = gmst + eq;

% ---------- 6. 极移矩阵 W ----------
% 精确旋转，等价于 Vallado 的 IAU-80 小角度矩阵
W = Ry(pm(1)) * Rx(pm(2));

% ---------- 7. 合成总矩阵: J2000 -> 地固系(ITRF) ----------
% 自右向左依次作用: 岁差 P -> 章动 N -> 自转 Rz(-GAST) -> 极移 W
C = W * Rz(-gast) * N * P;

% ---------- 8. 返回中间量 ----------
info = struct('T_TT', T, 'GMSTRad', gmst, 'GASTRad', gast, ...
    'dPsiRad', dpsi, 'dEpsilonRad', deps, 'P', P, 'N', N, 'W', W);
end