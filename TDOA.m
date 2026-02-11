%% ======================= TDOA (2D) 最少3基站：网格初值 + fminsearch 精化 + 双曲线可视化 =======================
clear; clc; close all;                                         % 清空变量/命令行/关闭图窗

c = 3e8;                                                        % 传播速度(米/秒)，无线电近似光速
sigma_dt = 0.5e-9;                                              % TDOA时间差测量噪声标准差(秒)，例如0.5ns

S = [0, 0;                                                       % 基站1坐标(参考站) [x1,y1]
     12, 0;                                                      % 基站2坐标 [x2,y2]
     0, 12];                                                     % 基站3坐标 [x3,y3]

p_true = [4.0; 5.0];                                             % 真实目标坐标 [x;y]

d_true = sqrt(sum((S - p_true.').^2, 2));                        % 真实距离 d_i = ||p - s_i||，得到3x1向量
dt21_true = (d_true(2) - d_true(1))/c;                           % 真实TDOA: t2 - t1 = (d2 - d1)/c
dt31_true = (d_true(3) - d_true(1))/c;                           % 真实TDOA: t3 - t1 = (d3 - d1)/c

dt21_meas = dt21_true + sigma_dt*randn;                          % 带噪TDOA测量(秒)
dt31_meas = dt31_true + sigma_dt*randn;                          % 带噪TDOA测量(秒)

Delta21 = c*dt21_meas;                                           % 距离差测量 Δ21 = d2 - d1 (米)
Delta31 = c*dt31_meas;                                           % 距离差测量 Δ31 = d3 - d1 (米)

s1 = S(1,:).';                                                   % 参考基站1坐标列向量 [x1;y1]
s2 = S(2,:).';                                                   % 基站2坐标列向量 [x2;y2]
s3 = S(3,:).';                                                   % 基站3坐标列向量 [x3;y3]

center = mean(S,1).';                                            % 用基站中心作为“合理区域”的中心(列向量)
Rmax = 200;                                                      % 给搜索一个最大半径(米)，防止优化发散到天边
w_pen = 1e-2;                                                    % 半径越界惩罚权重(越大越不容易跑远)

obj = @(p) ( (norm(p - s2) - norm(p - s1) - Delta21)^2 ...        % 目标函数：残差1平方
           + (norm(p - s3) - norm(p - s1) - Delta31)^2 ...        % 目标函数：残差2平方
           + w_pen * (max(0, norm(p - center) - Rmax))^2 );       % 软约束惩罚：超出Rmax才惩罚，防止发散

margin = 10;                                                     % 网格搜索边界外扩(米)
xmin = min(S(:,1)) - margin;                                     % 网格x最小值
xmax = max(S(:,1)) + margin;                                     % 网格x最大值
ymin = min(S(:,2)) - margin;                                     % 网格y最小值
ymax = max(S(:,2)) + margin;                                     % 网格y最大值

nx = 301;                                                        % 网格x方向点数(越大越准但越慢)
ny = 301;                                                        % 网格y方向点数(越大越准但越慢)

xg = linspace(xmin, xmax, nx);                                   % 生成x网格坐标
yg = linspace(ymin, ymax, ny);                                   % 生成y网格坐标
[X, Y] = meshgrid(xg, yg);                                       % 生成二维网格(ny行,nx列)

D1 = sqrt((X - S(1,1)).^2 + (Y - S(1,2)).^2);                    % 网格点到BS1距离
D2 = sqrt((X - S(2,1)).^2 + (Y - S(2,2)).^2);                    % 网格点到BS2距离
D3 = sqrt((X - S(3,1)).^2 + (Y - S(3,2)).^2);                    % 网格点到BS3距离

F21 = (D2 - D1) - Delta21;                                       % TDOA几何方程：||p-s2|| - ||p-s1|| - Δ21
F31 = (D3 - D1) - Delta31;                                       % TDOA几何方程：||p-s3|| - ||p-s1|| - Δ31

Jgrid = F21.^2 + F31.^2;                                         % 网格上的“无惩罚”残差平方(用于找初值)

[~, idx_min] = min(Jgrid(:));                                    % 找到残差平方最小的网格点索引
[y_idx, x_idx] = ind2sub(size(Jgrid), idx_min);                  % 将线性索引转换为(行,列)索引
p0 = [X(y_idx, x_idx); Y(y_idx, x_idx)];                         % 网格粗搜索得到的初值 p0=[x0;y0]

opts = optimset('Display','off', 'MaxIter', 2000, 'MaxFunEvals', 5000); % 设置fminsearch参数(不打印,给足迭代)
p_hat = fminsearch(obj, p0, opts);                               % 用Nelder-Mead从p0出发做非线性最小化

fprintf('TDOA: true = (%.3f, %.3f),  init = (%.3f, %.3f),  hat = (%.3f, %.3f)\n', ...
    p_true(1), p_true(2), p0(1), p0(2), p_hat(1), p_hat(2));     % 打印真实/初值/估计

figure; hold on; grid on; axis equal;                            % 新图窗、保持、网格、等比例坐标
title('TDOA (3 BS) - target location estimation');% 标题
xlabel('x (m)'); ylabel('y (m)');                                % 坐标轴标签

plot(S(:,1), S(:,2), 'ks', 'MarkerSize', 9, 'LineWidth', 1.5);   % 画基站：黑色方块
text(S(:,1)+0.2, S(:,2)+0.2, {'BS1(ref)','BS2','BS3'});          % 标注基站

plot(p_true(1), p_true(2), 'go', 'MarkerSize', 9, 'LineWidth', 2); % 真实目标：绿圈
plot(p0(1),     p0(2),     'bd', 'MarkerSize', 8, 'LineWidth', 2); % 初值：蓝菱形
plot(p_hat(1),  p_hat(2),  'rx', 'MarkerSize', 10, 'LineWidth', 2);% 估计目标：红叉

contour(X, Y, F21, [0 0], 'k--', 'LineWidth', 1.2);              % 画双曲线1：F21=0(虚线)
contour(X, Y, F31, [0 0], 'k--', 'LineWidth', 1.2);              % 画双曲线2：F31=0(虚线)

legend('Base Stations', 'True target', 'Grid init', 'Estimated target', 'TDOA locus (hyperbola)', ...
       'Location','best');                                       % 图例

xlim([xmin, xmax]); ylim([ymin, ymax]);                          % 显示范围与网格一致，便于观察
