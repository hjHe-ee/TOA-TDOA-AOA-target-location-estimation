%% ======================= AOA (2D) 最少2基站：两条方位射线交会/最小二乘 =======================
clear; clc; close all;                                       % 清空变量/命令行/关闭图窗

sigma_theta = deg2rad(1.0);                                   % 方位角测量噪声标准差(弧度)，例：1度

S = [0, 0;                                                     % 基站1坐标
     10, 2];                                                   % 基站2坐标
p_true = [4.5; 6.0];                                           % 真实目标坐标

theta_true = zeros(2,1);                                       % 存放真实方位角
theta_meas = zeros(2,1);                                       % 存放测得方位角(带噪声)

for i = 1:2                                                    % 对两个基站分别生成方位角测量
    dx = p_true(1) - S(i,1);                                   % 目标相对基站的x差
    dy = p_true(2) - S(i,2);                                   % 目标相对基站的y差
    theta_true(i) = atan2(dy, dx);                             % 真实方位角 atan2(y-y_i, x-x_i)
    theta_meas(i) = theta_true(i) + sigma_theta*randn;         % 测得方位角=真实+噪声
end

N = zeros(2,2);                                                % 线性约束矩阵 N*p = b (2条直线)
b = zeros(2,1);                                                % 右端 b

for i = 1:2                                                    % 对每条方位线构造法向约束
    ni = [-sin(theta_meas(i)); cos(theta_meas(i))];            % 射线方向u=[cos;sin]的法向n=[-sin;cos]
    N(i,:) = ni.';                                             % 第i行放n^T
    b(i)   = ni.' * S(i,:).';                                  % b_i = n^T*s_i，使 n^T(p - s_i)=0
end

p_hat = N\b;                                                   % 解两条直线交点(含噪时为最小二乘)

fprintf('AOA: true = (%.3f, %.3f),  hat = (%.3f, %.3f)\n', ...
    p_true(1), p_true(2), p_hat(1), p_hat(2));                 % 打印结果

figure; hold on; grid on; axis equal;                          % 新图窗、保持、网格、等比例
title('AOA (2 BS) - target location estimation');   % 标题
xlabel('x (m)'); ylabel('y (m)');                               % 坐标轴标签

plot(S(:,1), S(:,2), 'ks', 'MarkerSize', 9, 'LineWidth', 1.5);  % 画基站
text(S(:,1)+0.2, S(:,2)+0.2, {'BS1','BS2'});                    % 标注基站

plot(p_true(1), p_true(2), 'go', 'MarkerSize', 9, 'LineWidth', 2); % 真实目标
plot(p_hat(1),  p_hat(2),  'rx', 'MarkerSize', 10, 'LineWidth', 2);% 估计目标

L = 20;                                                        % 射线长度(用于画线段)
for i = 1:2                                                     % 画每个基站的方位射线(虚线)
    ux = cos(theta_meas(i));                                    % 射线方向x分量
    uy = sin(theta_meas(i));                                    % 射线方向y分量
    x_line = [S(i,1), S(i,1) + L*ux];                            % 射线端点x坐标(从基站出发)
    y_line = [S(i,2), S(i,2) + L*uy];                            % 射线端点y坐标(从基站出发)
    plot(x_line, y_line, 'k--', 'LineWidth', 1.0);               % 虚线画射线
end

legend('Base Stations','True target','Estimated target','Bearing rays', ...
    'Location','best');                                         % 图例

xlim([-4, 18]); ylim([-4, 18]);                                  % 设置显示范围
