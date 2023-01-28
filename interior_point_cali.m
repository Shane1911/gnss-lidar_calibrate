close all
clear
clc
format long
load data/relative.mat
% ------------------------data_deal-mydata--------------------
% date: x, y, z, yaw, pitch, roll
%% Creat V and W;
[num, d] = size(Lidar_pose);
MA = zeros(4, 4 * num);
MB = zeros(4, 4 * num);
MA_d = zeros(4, 4 * num);
MB_d = zeros(4, 4 * num);
V = zeros(4, 4 * num);
W = zeros(4, 4 * num);
Q = zeros(8, 8);
for i = 1:num
    qr_a = eul2quat(Ins_pose(i, 4:6));
    qr_b = eul2quat(Lidar_pose(i, 4:6));
    % get V
    MA(:, 4 * i - 3:4 * i) = Trans_Ml(qr_a);
    MB(:, 4 * i - 3:4 * i) = Trans_Mr(qr_b);
    V(:, 4 * i - 3:4 * i) = MA(:, 4 * i - 3:4 * i) - MB(:, 4 * i - 3:4 * i);
    % get W
    qt_a = [0 Ins_pose(i, 1:3)];
    qt_b = [0 Lidar_pose(i, 1:3)];
    MA_d(:, 4 * i - 3:4 * i) = Trans_Ml(0.5 * quatmultiply(qt_a, qr_a));
    MB_d (:, 4 * i - 3:4 * i) = Trans_Mr(0.5 * quatmultiply(qt_b, qr_b));
    W(:, 4 * i - 3:4 * i) = MA_d(:, 4 * i - 3:4 * i) - MB_d(:, 4 * i - 3:4 * i);
    H1 = zeros(8, 8);
    H1(1:4, 1:4) = V(:, 4 * i - 3:4 * i);
    H1(5:8, 5:8) = V(:, 4 * i - 3:4 * i);
    H1(5:8, 1:4) = W(:, 4 * i - 3:4 * i);
    Q = Q + H1' * H1;
end
Q = 0.5 * (Q' + Q);
rng default

% For reproducibility
H{1} = [eye(4) zeros(4, 4); zeros(4, 4) zeros(4, 4)];
H2 = [zeros(4, 4) eye(4); zeros(4, 4) zeros(4, 4)];
H{2} = 0.5 * (H2 + H2');
d1{1} = -1;
d1{2} = 0;
% options = optimoptions(@fmincon, 'Algorithm', 'interior-point');
options = optimoptions(@fmincon, 'Algorithm', 'interior-point', ...
'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
    'HessianFcn', @(x, lambda)quadhess(x, lambda, Q, H));
fun = @(x)quadobj(x, Q);
nonlconstr = @(x)quadconstr(x, H, d1);
% q_r = eul2quat([0.001 0.001 0.001], 'ZYX');
% q_t = 0.5 * quatmultiply([0 -8 0.1 0.1], q_r);
% q_dul = [q_r q_t];
% x0 = q_dul; % Column vector
% x0 = rand(8, 1);
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
% lb = [-1 -1 -1 -1 -1 -1 -1 -1];
% ub = [1 1 1 1 1 1 1 1];
% lb = [-1, -1, -1, -pi, -pi, -pi];
% ub = [1, 1, 1, pi, pi, pi];
% x0 = randn(1, 8);
q_r = eul2quat([0.012 0.015 -0.11], 'ZYX');
q_t = 0.5 * quatmultiply([0 9.5 -0.6 0.35], q_r);
q_dul = [q_r q_t];
x0 = q_dul; % Column vector
% x0 =[0.697185788758996	-0.0182430780369156	0.0211713399555684	0.716345545420194	-0.251197377399799	1.90320139090422	1.78641034061645	0.240150572326555];
[x1, fval, eflag, output, lambda] = fmincon(fun, x0', ...
A, b, Aeq, beq, lb, ub, nonlconstr, options);
lambda.eqnonlin
x = zeros(1, 6);
x_r = x1(1:4)';
x(4:6) = quat2eul(x1(1:4)', 'ZYX');
x_t = 2 * quatmultiply(x1(5:8)', [x_r(1) -x_r(2:4)]);
x(1:3) = x_t(2:4);
% x = x1;
T12 = eul2tform(x(4:6), 'ZYX');
T12(1:3, 4) = x(1:3);
% T12 = inv(T12);
fprintf("T12 = \n")
disp(T12)
T2eul = [T12(1:3, 4)', rotm2eul(T12(1:3, 1:3), 'ZYX')];
save result/res.mat T2eul

function hess = quadhess(x, lambda, Q, H)
    hess = Q;
    jj = length(H); % jj is the number of inequality constraints

    for i = 1:jj
        hess = hess + lambda.eqnonlin(i) * H{i};
    end

end

function [y, grady] = quadobj(x, Q)
    % q_r = eul2quat(x(1, 4:6), 'ZYX');
    % q_t = 0.5 * quatmultiply([0 x(1:3)], q_r);
    % q_dul = [q_r q_t];
    y = 1/2 * x' * Q * x;

    if nargout > 1
        grady = Q * x;
    end

end

function [yeq, y, gradyeq, grady] = quadconstr(x, H, d)
    jj = length(H); % jj is the number of inequality constraints
    y = zeros(1, jj);

    for i = 1:jj
        y(i) = x' * H{i} * x + d{i};
    end

    yeq = [];

    if nargout > 2
        grady = zeros(length(x), jj);

        for i = 1:jj
            grady(:, i) = H{i} * x;
        end

    end

    gradyeq = [];
end

function M = Trans_Ml(quat_r)
    M = [quat_r(1) -quat_r(2) -quat_r(3) -quat_r(4);
        quat_r(2) quat_r(1) -quat_r(4) quat_r(3);
        quat_r(3) quat_r(4) quat_r(1) -quat_r(2);
        quat_r(4) -quat_r(3) quat_r(2) quat_r(1)];
end

function M = Trans_Mr(quat_r)
    M = [quat_r(1) -quat_r(2) -quat_r(3) -quat_r(4);
        quat_r(2) quat_r(1) quat_r(4) -quat_r(3);
        quat_r(3) -quat_r(4) quat_r(1) quat_r(2);
        quat_r(4) quat_r(3) -quat_r(2) quat_r(1)];
end

