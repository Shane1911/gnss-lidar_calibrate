close all
clear
clc
% gnss_time = load("data/time_gps.txt");
gnss_time = load("data_back_cali/localition_data/time_gps.txt");
% x,y,z,ve,vn,vu,qx,qy,qz,qw
% gnss_pose = load("data/gps-0.txt");
gnss_pose = load("data_back_cali/localition_data/gps_data.txt");
% cloud_time = load("data/time_cloud.txt");
cloud_time = load("data_back_cali/localition_data/time_cloud.txt");
time0 = gnss_time(1);
n = find(cloud_time  > time0);
cloud_time = cloud_time(n);
[time_lidar_filter,m,~]= unique(cloud_time,'rows');
% k-:the index of gps which is nearest one away from the lidar pose.
k = dsearchn(gnss_time,time_lidar_filter);

% tg = gnss_time(1:20);
% dt = abs(tg - cloud_time(1));
% list_id -> get pose:x,y,z,qw,qx,qy,qz
list_id = [1, 2, 3, 10, 7, 8, 9];
pose_gps = gnss_pose(:,list_id);
% way-1:nearst_pose(no slerp)
nearest_gps = pose_gps(k,:);
figure
plot(nearest_gps(:, 1), nearest_gps(:, 2),  'rs-', 'LineWidth', 1)
% hold on 
% plot(ex(:,1),ex(:,2),'b.-');
xlabel('X / m')
ylabel('Y / m')
% zlabel('Z / m')
title('最近点')
legend('Gps-lidar odometry')

% way-2:slerp
nums = size(cloud_time);
slerp_pose = [];
for i = 1:nums
    dt = cloud_time(i) - gnss_time(k(i));
    if dt > 0
        id_left = k(i);
        id_right = k(i)+1;
    else
        id_right = k(i);
        id_left = k(i)-1;
    end
    lambda = (cloud_time(i) - gnss_time(id_left))/(gnss_time(id_right)-gnss_time(id_left));
    slerp_pose(i,:) = slerp_gps(id_left,id_right,lambda,pose_gps);
end
save data/nearest_gps.txt -ascii nearest_gps;
save data/slerp_gps.txt -ascii slerp_pose;
figure
plot(slerp_pose(:, 1), slerp_pose(:, 2),  'rs-', 'LineWidth', 2)
xlabel('X / m')
ylabel('Y / m')
% zlabel('Z / m')
title('插值点')
legend('Gps-lidar odometry')

% compare
d_pose = slerp_pose(:,1:3)-nearest_gps(:,1:3);
d_eul = quat2eul(slerp_pose(:,4:7)) - quat2eul(nearest_gps(:,4:7));
figure
subplot(2,3,1)
plot(d_pose(:,1),  'rs-')
title('dx')
subplot(2,3,2)
plot(d_pose(:,2),  'rs-')
title('dy')
subplot(2,3,3)
plot(d_pose(:,3),  'rs-')
title('dz')
subplot(2,3,4)
plot(d_eul(:,1),  'bs-')
title('yaw')
subplot(2,3,5)
plot(d_eul(:,2),  'bs-')
title('pitch')
subplot(2,3,6)
plot(d_eul(:,3),  'bs-')
title('roll')

function t_slerp_pose = slerp_gps(t_l,t_r,l,pose_gps)
    % pose slerp
    pose = (1-l).*pose_gps(t_l,1:3) + l.*pose_gps(t_r,1:3);
    % q slerp
    quat = slerpQuaternions(pose_gps(t_l,4:7),pose_gps(t_r,4:7),l);
    t_slerp_pose = [pose quat];
end

function [ Qslerped ] = slerpQuaternions ( Q1, Q2, W )
    q1 = quaternion(Q1);
    q2 = quaternion(Q2);
    q = slerp(q1,q2,W);
    angle_eul = euler(q,"ZYX","frame");
    Qslerped = eul2quat(angle_eul);
end
