%% 输入六根数及 UTC 时刻，计算该时刻的 J2000 位置和速度
% 直接修改下面的输入区，然后运行本脚本。
% 第六根数为平近点角 M；半长轴为 m，角度为 deg。
% 输出 rJ2000_m、vJ2000_mps 均为 3×1 列向量。
clc; clear; 
%% 1. 输入区
%          a(m)       e    i(deg)   RAAN(deg)  w(deg)  M(deg)
coe_m = [6922975      0.001    39.123  354.453    95.592  63.363];

%             年    月  日  时  分  秒（保留三位小数），UTC
epoch_utc = [2026   10  1   4  55  22.190];

% STK Coordinate Epoch：通常与 Orbit Epoch 相同。
% 若 STK 中设置了不同的 Coordinate Epoch，请单独修改此项。
coordinate_epoch_utc = epoch_utc;

% 默认读取本脚本所在目录中的 EOP 文件，也可改为完整路径。
script_dir = fileparts(mfilename('fullpath'));
eop_file = fullfile(script_dir, 'EOP-Last5Years.txt');

%% 2. 计算
validateattributes(coe_m, {'numeric'}, {'real','finite','vector','numel',6});
coe_km = double(coe_m(:).');
coe_km(1) = coe_km(1) / 1000;  % 核心函数内部使用 km、km/s

[rJ_km, vJ_kmps, info] = AlignmentAtEpoch_to_J2000( ...
    coe_km, epoch_utc, eop_file, ...
    'CoordinateEpoch', coordinate_epoch_utc);

rJ2000_m = rJ_km * 1000;
vJ2000_mps = vJ_kmps * 1000;
rvJ2000 = [rJ2000_m; vJ2000_mps];  % [x;y;z;vx;vy;vz]，m、m/s

%% 3. 输出
t = info.OrbitEpochUTC;
fprintf('UTC: %04d-%02d-%02d %02d:%02d:%06.3f\n', t);
fprintf('J2000 r (m)   : [% .6f  % .6f  % .6f]\n', rJ2000_m);
fprintf('J2000 v (m/s) : [% .6f  % .6f  % .6f]\n', vJ2000_mps);
fprintf('EOP: %s\n', info.EOP.File);
if info.EOP.IsPredicted
    fprintf('该时刻使用预测 EOP 数据。\n');
end
% stk -2802232.284    -6117449.779    1623932.926    5056.839646    -3497.319197    -4447.825060