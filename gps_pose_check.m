clear
clc
format long
imu_time = load("data_back_cali/localition_data/time_imu.txt");
imu_data = load("data_back_cali/localition_data/imu_data.txt");
lidar_time = load("data_back_cali/localition_data/time_cloud.txt");
% load
data_lidar = load('data_back_cali/lidar_odometry.txt');
% x,y,z,qw,qx,qy,qz
list_id = [5,6,7,1,2,3,4];
Lidar_pose = data_lidar(:,list_id);
% LL2UTM
[gps_utm] = 1e3*lat_lon2utm(imu_data(:,1),imu_data(:,2));
gps_utm =[gps_utm imu_data(:,3)];
gnss_pose = load("data_back_cali/localition_data/gps_data.txt");
gnss_pose = gnss_pose(1:end,:);
% get eulangles
eul_imu = imu_data(:,[7 9 8])*pi/180;
eul_imu(:,1) =  pi/2 - eul_imu(:,1);
q1 =quatnormalize(eul2quat(eul_imu));
eul_imu = quat2eul(q1);
% check dulangles
v_enu = imu_data(:,4:6);
theta_enu = [];
for i = 1:length(v_enu(:,1))
    theta = atan(v_enu(i,2)./v_enu(i,1))*180/pi;
    if (v_enu(i,1)<0 && v_enu(i,2)<0)
        theta = theta - 180;
    end
    if (v_enu(i,1)<0 && v_enu(i,2)>0)
        theta = theta + 180;
    end
    theta_enu = [theta_enu;theta];
end
% view the trajectory(right-east,up-north)
figure
subplot(1,2,1)
plot(theta_enu,'r.-',LineWidth=1)
hold on 
plot(eul_imu(:,1)*180/pi,'b.-',LineWidth=1)
grid on 
legend('eul from venu','eul from pose')
subplot(1,2,2)
% gps_utm(:,1:2) = gps_utm(:,1:2) - gps_utm(1,1:2);
plot(gps_utm(:, 1), gps_utm(:, 2),  'r.-', LineWidth=1)
hold on 
% gnss_pose(:,1:3) = gnss_pose(:,1:3) - gnss_pose(1,1:3);
plot(gnss_pose(:, 1), gnss_pose(:, 2),  'b.-', LineWidth=1)
xlabel('X / m')
ylabel('Y / m')
% zlabel('Z / m')
title('最近点')
legend('ll2utm','Gps-lidar odometry')

% get the relative pose
% time 
time0 = imu_time(1);
n = lidar_time  > time0;
lidar_time = lidar_time(n);
[time_lidar_filter,~,~]= unique(lidar_time,'rows');
% k-:the index of gps which is nearest one away from the lidar pose.
k = dsearchn(imu_time,time_lidar_filter);
% knum = length(k);
knum = length(find(k<max(length(gps_utm),length(eul_imu))));
gps_utm = gps_utm(k(1:knum),:);
eul_imu = eul_imu(k(1:knum),:);


% --gnss/lidar:相对位置
R0_1 = quat2rotm(Lidar_pose(1, 4:7)); % qw qx qy qz
t0_1 = Lidar_pose(1, 1:3);
[m, ~] = size(Lidar_pose);
for i = 1:m
    Lidar_pose(i, 1:3) = R0_1 \ (Lidar_pose(i, 1:3)' - t0_1');
    R = R0_1 \ quat2rotm(Lidar_pose(i, 4:7)); % qw qx qy qz
    eul = rotm2eul(R, 'ZYX'); % rad
    Lidar_pose(i, 4:6) = eul; % ZYX
end
Lidar_pose(:, 7) = [];

R0_2 = eul2rotm(eul_imu(1, :)); % qw qx qy qz
t0_2 = gps_utm(1, 1:3);
[n, ~] = size(gps_utm(:,1));
Ins_pose = [];
for i = 1:n
    Ins_pose(i, 1:3) = R0_2 \ (gps_utm(i, 1:3)' - t0_2');
    R = R0_2 \ eul2rotm(eul_imu(i, :)); % qw qx qy qz
    eul = rotm2eul(R, 'ZYX'); % rad
    Ins_pose(i, 4:6) = eul; % ZYX
end
figure
plot(Ins_pose(:,1),Ins_pose(:,2),'r.-',LineWidth=1)
hold on 
plot(Lidar_pose(:,1),Lidar_pose(:,2),'b.-',LineWidth=1)
grid on 
legend('gnss pose','lidar')
save data/relative.mat Ins_pose Lidar_pose

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
