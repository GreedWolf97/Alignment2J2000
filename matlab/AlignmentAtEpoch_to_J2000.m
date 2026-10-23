function [rJ, vJ, info] = AlignmentAtEpoch_to_J2000(coe, epoch, EOP_FILE, varargin)
% coe=[a(km),e,i(deg),RAAN(deg),argPerigee(deg),M(deg)]. UTC epoch.
% Returns state at orbit epoch, not a propagated trajectory.
% CoordinateEpoch defaults to orbit epoch. No Aerospace Toolbox required.
% Options: Mu (398600.4415 km^3/s^2), UsePolarMotion (true),
% UseNutationCorrections (false for STK J2000). rJ in km, vJ in km/s.

% 若未指定 EOP 文件，则默认使用与本函数同目录下的 EOP-Last5Years.txt
if nargin < 3 || isempty(EOP_FILE)
    EOP_FILE = fullfile(fileparts(mfilename('fullpath')), 'EOP-Last5Years.txt');
end

% 校验轨道六根数：必须是 6 个实数有限值
validateattributes(coe, {'numeric'}, {'real', 'finite', 'vector', 'numel', 6});

% 校验并规范化历元时间（UTC [Y M D h m s]）
epoch = valid_epoch(epoch);

% 解析可选参数
p = inputParser;
addParameter(p, 'CoordinateEpoch', epoch);   % 坐标历元，默认等于轨道历元
% STK 12.2 Earth.cb Gm = 3.986004415e14 m^3/s^2.
addParameter(p, 'Mu', 398600.4415, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0); % 地球引力常数
addParameter(p, 'UsePolarMotion', true, @(x) isscalar(x) && any(x == [0 1]));               % 是否使用极移
% STK J2000 uses uncorrected FK5 nutation. EOP dPsi/dEpsilon belong
% to its J2000/ICRF relation; adding them here changes this frame.
addParameter(p, 'UseNutationCorrections', false, @(x) isscalar(x) && any(x == [0 1]));      % 是否加入章动改正
parse(p, varargin{:});
o = p.Results;

% 校验坐标历元
ce = valid_epoch(o.CoordinateEpoch);

% 提取轨道六根数并转为行向量
c = double(coe(:).');
a = c(1);   % 半长轴 (km)
e = c(2);   % 偏心率

% 根数合法性检查
if a <= 0 || e < 0 || e >= 1 || c(3) < 0 || c(3) > 180
    error('Alignment:InvalidElements', 'Requires a>0, 0<=e<1, 0<=i<=180.');
end

% 读取 EOP 数据（极移 xp/yp、UT1-UTC、dAT、章动 dPsi/dEpsilon）
q = read_alignment_eop(EOP_FILE, ce);
as2r = pi / 648000;   % 角秒 -> 弧度

% 极移（弧度），乘以开关；章动改正同理（弧度）
pm = [q.xp q.yp] * as2r * double(o.UsePolarMotion);
dn = [q.dPsi q.dEpsilon] * as2r * double(o.UseNutationCorrections);

% 计算 J2000 -> 地固系 的旋转矩阵 C，以及帧相关的辅助信息
[C, frameInfo] = fk5_j2000_to_fixed(ce, q.deltaAT, q.deltaUT1, pm, dn);

% 地固系 -> J2000 的旋转矩阵（转置）
A = C.';

% 平近点角归一化到 [-pi, pi)
M = mod(c(6) * pi / 180 + pi, 2 * pi) - pi;
lo = -pi;
hi = pi;

% Monotonic Kepler equation, robust for high eccentricity.
% 二分法求解开普勒方程 M = E - e*sin(E)，对高偏心率也稳定
for k = 1:60
    E = (lo + hi) / 2;
    if E - e * sin(E) > M
        hi = E;
    else
        lo = E;
    end
end
E = (lo + hi) / 2;   % 收敛后的偏近点角

% 由偏近点角计算近焦点坐标系(PQW)下的位置、速度
b = sqrt(1 - e^2);
rp = a * [cos(E) - e; b * sin(E); 0];                                  % 位置 (km)
vp = sqrt(o.Mu / a) / (1 - e * cos(E)) * [-sin(E); b * cos(E); 0];     % 速度 (km/s)

% 定义绕 Z 轴、X 轴的旋转矩阵函数
Rz = @(x) [cos(x) -sin(x) 0; sin(x) cos(x) 0; 0 0 1];
Rx = @(x) [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];

% 构造 PQW -> 地固系 的旋转矩阵：Rz(RAAN) * Rx(i) * Rz(argPerigee)
Q = Rz(c(4) * pi / 180) * Rx(c(3) * pi / 180) * Rz(c(5) * pi / 180);

% 将位置、速度转到地固系（Alignment 坐标系）
rA = Q * rp;
vA = Q * vp;

% AlignmentAtEpoch is frozen: no Earth angular velocity cross r term.
% 说明：本函数冻结在历元时刻，不考虑地球自转引起的附加速度项
% 由地固系转到 J2000
rJ = A * rA;
vJ = A * vA;

% 汇总所有中间量到 info 结构体，便于事后检查与调试
info = struct('OrbitEpochUTC', epoch, 'CoordinateEpochUTC', ce, 'Mu', o.Mu, ...
    'Model', 'IAU-76/FK5', 'EOP', q, 'rAlignment', rA, 'vAlignment', vA, ...
    'C_J2000_to_Fixed', C, 'C_Fixed_to_J2000', A, 'Options', o, 'Frame', frameInfo);
end

function t = valid_epoch(t)
% 校验历元格式：6 个实数有限值，且符合 [Y M D h m s] 的合理范围
validateattributes(t, {'numeric'}, {'real', 'finite', 'vector', 'numel', 6});
t = double(t(:).');

% 年/月/日/时/分必须为整数；秒可以为小数
% 月 1..12，日 >=1，时 0..23，分 0..59，秒 0..<60
if any(t(1:5) ~= fix(t(1:5))) || t(1) < 1 || t(2) < 1 || t(2) > 12 || t(3) < 1 || t(4) < 0 || t(4) >= 24 || t(5) < 0 || t(5) >= 60 || t(6) < 0 || t(6) >= 60
    error('Alignment:InvalidEpoch', 'Expected UTC [Y M D h m s], 0<=s<60.');
end

% 使用 datetime 校验日期是否为合法日历日期（如非闰年 2 月 30 日应报错）
d = datetime(t(1), t(2), t(3));
if year(d) ~= t(1) || month(d) ~= t(2) || day(d) ~= t(3)
    error('Alignment:InvalidEpoch', 'Invalid calendar date.');
end
end