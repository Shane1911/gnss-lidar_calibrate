close all
clear
clc
format long
load data/relative.mat
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
%way-1 SDP_TRM
Q0 = 2*Q;
Eqcon(1).Q = zeros(8, 8); %H1
Eqcon(1).Q(1:4, 1:4) = eye(4);
Eqcon(1).c = -1;
%weight_factor
limit_mu = 0;
H2 = [limit_mu*eye(4) eye(4); zeros(4, 4) zeros(4, 4)];
H_2 = 0.5*(H2 + H2');
Eqcon(2).Q = H_2; %H1
Eqcon(2).c = -1*limit_mu;
Incon = [];
output = irma(Q0, Incon, Eqcon);
x1 = output.x;
% angle = quat2eul(x1(1:4)', 'ZYX') * 180 / pi
pose_x = zeros(1, 6);
x_r = x1(1:4);
pose_x(4:6) = quat2eul(x1(1:4)', 'ZYX');
x_t = 2 * quatmultiply(x1(5:8)', [x_r(1); -x_r(2:4)]');
pose_x(1:3) = x_t(2:4);
T2eul = pose_x;
save result/res.mat T2eul
disp(x_t)
disp(pose_x(4:6)*180/pi)

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
