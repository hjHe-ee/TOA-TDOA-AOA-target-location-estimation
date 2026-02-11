%% ======================= TOA (2D) 最少3基站：圆交会/最小二乘 =======================
clear; clc; close all;                                  % 清空变量/命令行/关闭图窗

c = 3e8;                                                 % 传播速度(无线电近似光速)
sigma_r = 0.30;                                          % 距离测量噪声标准差(米)

S = [0, 0;                                                % 基站1坐标 [x1,y1]
     10, 0;                                               % 基站2坐标 [x2,y2]
     0, 10];                                              % 基站3坐标 [x3,y3]
p_true = [3.2; 4.6];                                      % 真实目标坐标 [x;y]

d_true = sqrt(sum((S - p_true.').^2, 2));                 % 真实距离向量 d_i = ||p - s_i||
r_meas = d_true + sigma_r*randn(size(d_true));            % 测得距离(TOA换算得到的距离)，加高斯噪声

x1 = S(1,1); y1 = S(1,2);                                 % 取基站1坐标，便于写公式
r1 = r_meas(1);                                           % 取基站1测得距离

A = zeros(2,2);                                           % 构造线性方程 A*[x;y]=b (2条方程,2未知)
b = zeros(2,1);                                           % 构造线性方程右端 b

for k = 2:3                                               % 用基站2、3分别与基站1作差，得到线性方程
    xi = S(k,1); yi = S(k,2);                             % 第k个基站坐标
    ri = r_meas(k);                                       % 第k个基站测得距离
    A(k-1,:) = [2*(xi - x1), 2*(yi - y1)];                % 2(xi-x1)x + 2(yi-y1)y = ...
    b(k-1)   = r1^2 - ri^2 + xi^2 - x1^2 + yi^2 - y1^2;   % ... = r1^2-ri^2 + (xi^2-x1^2)+(yi^2-y1^2)
end

p_hat = A\b;                                              % 最小二乘/直接解(2x2时为精确解)，得到估计位置

fprintf('TOA: true = (%.3f, %.3f),  hat = (%.3f, %.3f)\n', ...
    p_true(1), p_true(2), p_hat(1), p_hat(2));            % 打印结果

figure; hold on; grid on; axis equal;                     % 新图窗、保持、网格、等比例坐标
title('TOA (3 BS) - target location estimation');  % 图标题
xlabel('x (m)'); ylabel('y (m)');                         % 坐标轴标签

plot(S(:,1), S(:,2), 'ks', 'MarkerSize', 9, 'LineWidth', 1.5); % 画基站：黑色方块
text(S(:,1)+0.2, S(:,2)+0.2, {'BS1','BS2','BS3'});        % 标注基站编号

plot(p_true(1), p_true(2), 'go', 'MarkerSize', 9, 'LineWidth', 2); % 真实目标：绿圈
plot(p_hat(1),  p_hat(2),  'rx', 'MarkerSize', 10, 'LineWidth', 2);% 估计目标：红叉

theta = linspace(0, 2*pi, 400);                           % 圆参数角度
for i = 1:3                                               % 依次画每个基站的“测量圆”
    cx = S(i,1) + r_meas(i)*cos(theta);                   % 圆的x坐标：xi + ri*cos(t)
    cy = S(i,2) + r_meas(i)*sin(theta);                   % 圆的y坐标：yi + ri*sin(t)
    plot(cx, cy, 'k--', 'LineWidth', 1.0);                % 虚线圆（测量几何）
end

legend('Base Stations','True target','Estimated target','Measured circles', ...
    'Location','best');                                   % 图例

xlim([-2, 14]); ylim([-2, 14]);                           % 设置显示范围，便于观察
