close all
clear
clc
%====================Load original data===============
format long
data_gnss = load('data/slerp_gps.txt');
data_lidar = load('data_back_cali/lidar_odometry.txt');
num = min(length(data_gnss),length(data_lidar));
list_id = [5,6,7,1,2,3,4];
slerp_pose = data_gnss(1:num,:);
Lidar_pose = data_lidar(1:num,list_id);
save data/origional.mat Lidar_pose slerp_pose

% --gnss/lidar:相对位置
% slip sample for correct the error.
len = length(Lidar_pose(:,1));
% the length of slip.
delta_num = 50;
num_split = floor(len/delta_num);
indx = [];
Ins_pose1 = [];
Lidar_pose1 = [];
for j = 1:num_split
    indx = ((1+delta_num*(j-1)):(delta_num*j));
    [Ins_posej,Lidar_posej] =getrelativepose(Lidar_pose(indx,:),slerp_pose(indx,:));
    Ins_pose1 = [Ins_pose1;Ins_posej];
    Lidar_pose1 = [Lidar_pose1;Lidar_posej];
end
[Ins_posej,Lidar_posej] =getrelativepose(Lidar_pose(indx(end):end,:),slerp_pose(indx(end):end,:));
Ins_pose1 = [Ins_pose1;Ins_posej];
Lidar_pose1 = [Lidar_pose1;Lidar_posej];

%get scan to scan pose
Tl_pre = [quat2rotm(Lidar_pose(1, 4:7)) Lidar_pose(1, 1:3)';0 0 0 1];
[m, ~] = size(Lidar_pose);
re_lidar_pose = size(m,7);
delta_id = (1:m);
for k = 1:length(delta_id)
    i = delta_id(k);
    Tl_cur = [quat2rotm(Lidar_pose(i, 4:7)) Lidar_pose(i, 1:3)';0 0 0 1];
    dt = Tl_pre \ Tl_cur;
    re_lidar_pose(k, 4:6) = rotm2eul(dt(1:3,1:3), 'ZYX'); % ZYX
    re_lidar_pose(k, 1:3) = dt(1:3,4)'; % ZYX
    Tl_pre = Tl_cur;
end

Tg_pre = [quat2rotm(slerp_pose(1, 4:7)) slerp_pose(1, 1:3)';0 0 0 1];
[n, ~] = size(slerp_pose);
re_ins_pose = size(n,7);
delta_id = (1:n);
for k = 1:length(delta_id)
    i = delta_id(k);
    Tg_cur = [quat2rotm(slerp_pose(i, 4:7)) slerp_pose(i, 1:3)';0 0 0 1];
    dt = Tg_pre \ Tg_cur;
    re_ins_pose(k, 4:6) = rotm2eul(dt(1:3,1:3), 'ZYX'); % ZYX
    re_ins_pose(k, 1:3) = dt(1:3,4)'; % ZYX
    Tg_pre = Tg_cur;
end
[num_front, ~] = size(Lidar_pose1);
select_front_id = (1:1:num_front);
select_back_id = (1:600);
Lidar_pose_A = [Lidar_pose1;re_lidar_pose(select_back_id,:)];
Ins_pose_B = [Ins_pose1;re_ins_pose(select_back_id,:)];
Ins_pose = Ins_pose_B;
Lidar_pose = Lidar_pose_A;
save data/relative.mat Ins_pose Lidar_pose

function [ins_t,lidar_t] =getrelativepose(Lidar_pose1,slerp_pose1)
    R0_1 = quat2rotm(Lidar_pose1(1, 4:7)); % qw qx qy qz
    t0_1 = Lidar_pose1(1, 1:3);
    [m, ~] = size(Lidar_pose1);
    lidar_t = [];
    for i = 1:m
        lidar_t(i, 1:3) = R0_1 \ (Lidar_pose1(i, 1:3)' - t0_1');
        R = R0_1 \ quat2rotm(Lidar_pose1(i, 4:7)); % qw qx qy qz
        eul = rotm2eul(R, 'ZYX'); % rad
        lidar_t(i, 4:6) = eul; % ZYX
    end

    R0_2 = quat2rotm(slerp_pose1(1, 4:7)); % qw qx qy qz
    t0_2 = slerp_pose1(1, 1:3);
    [n, ~] = size(slerp_pose1);
    ins_t = [];
    for i = 1:n
        ins_t(i, 1:3) = R0_2 \ (slerp_pose1(i, 1:3)' - t0_2');
        R = R0_2 \ quat2rotm(slerp_pose1(i, 4:7)); % qw qx qy qz
        eul = rotm2eul(R, 'ZYX'); % rad
        ins_t(i, 4:6) = eul; % ZYX
    end

end
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

function[shuchu]=lat_lon2utm(lat_shuru,lon_shuru)
    %地理经纬度坐标转换为UTM坐标
    %size_shuzu=size(lat_shuru);
    for i=1:length(lat_shuru)
        %输入经纬度
        lat=lat_shuru(i);
        lon=lon_shuru(i);
        %%unit：km
        a=6378.137; 
        e=0.0818192;
        k0=0.9996;
        E0=500;
        N0=0;
        Zonenum=fix(lon/6)+31;
        lamda0=(Zonenum-1)*6-180+3;%degree
        lamda0=lamda0*pi/180;%radian
        phi=lat*pi/180;
        lamda=lon*pi/180;%radian
        v=1/sqrt(1-e^2*(sin(phi)^2));
        A=(lamda-lamda0)*cos(phi);
        T=tan(phi)*tan(phi);
        C=e^2*cos(phi)*cos(phi)/(1-e^2);
        s=(1-e^2/4-3*e^4/64-5*e^6/256)*phi-(3*e^2/8+3*e^4/32+45*e^6/1024)*sin(2*phi)+...
            (15*e^4/256+45*e^6/1024)*sin(4*phi)-35*e^6/3072*sin(6*phi);
        UTME=E0+k0*a*v*(A+(1-T+C)*A^3/6+(5-18*T+T^2)*A^5/120);
        UTMN=N0+k0*a*(s+v*tan(phi)*(A^2/2+(5-T+9*C+4*C^2)*A^4/24+(61-58*T+T^2)*A^6/720));
        shuchu(i,1)=UTME;
        shuchu(i,2)=UTMN;
    end
end
