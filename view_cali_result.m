close all
clear
clc
format long
load result/res.mat
disp('==================================Cali R and t====================================')
disp(T2eul)
load data/origional.mat

T12 = eye(4);
% T2eul(2) = 0;
T12(1:3,1:3) = eul2rotm(T2eul(4:6));
T12(1:3,4) = T2eul(1:3);
start_id = 300;
T_g0 = get_pose(slerp_pose(start_id,:));
end_id = start_id + 200;
T_l0 = get_pose(Lidar_pose(start_id,:));
Lidar2gnss = [];
for j = start_id:end_id
    pose_l = eye(4);
    pose_l(1:3,4) = T_l0(1:3,1:3) \ (Lidar_pose(j,1:3)' - T_l0(1:3,4));
    pose_l(1:3,1:3) = T_l0(1:3,1:3) \ quat2rotm(Lidar_pose(j, 4 : 7)); % qw qx qy qz
%     pose_l = T_l0 \ get_pose(Lidar_pose(j,:));
    T_l2g = T_g0*(T12 *pose_l /T12);
    q = rotm2quat(T_l2g(1:3,1:3));
    t = T_l2g(1:3,4);
    Lidar2gnss = [Lidar2gnss;t' q];
end
Lidar2gnss = Lidar2gnss - slerp_pose(1,:);
slerp_pose = slerp_pose - slerp_pose(1,:);
figure
grid on
axis equal
plot(Lidar2gnss(:, 1), Lidar2gnss(:, 2),'ks-')
hold on
plot(slerp_pose(:, 1), slerp_pose(:, 2), 'r.-')
xlabel('X / m')
ylabel('Y / m')
zlabel('Z / m')
zlim([-20 30])
title('标定前')
legend('LiDAR Original', 'INS')

function pose = get_pose(q)
    T12 = eye(4);
    T12(1:3,1:3) = quat2rotm(q(4:7));
    T12(1:3,4) = q(1:3);
    pose = T12;
end
