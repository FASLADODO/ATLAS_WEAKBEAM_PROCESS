%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%Copyright(c) 2021, Optoelectric Detection and Laser Remote Sensing Laboratory, WuHanUniversity, CHINA
% Project: Space-boren LiDAR Simulation
% Team leader: Prof. Song Li, A/Prof. Yue Ma
% Developer: Zhiyu Zhang, Xinyuan Liu
% Contact Info:Zhiyuzhang@whu.edu.cn, mayue19860103@163.com; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear

%% Load data and data pre-processing
load('data_20190211035235.mat');%load geolocation photon data of strong beam
TBLHC1 = TBLHC_2th_strong;
PE_distance = TBLHC1(:,1) - TBLHC1(1,1);
PE_deltatime = TBLHC1(:,6);
PE_height = TBLHC1(:,2);
conf  = TBLHC1(:,5);
TBLHC_PE = [PE_distance,PE_height,PE_deltatime,conf];%save the processed data

load('data_20190211035235_dbscan.mat');%load DBSCAN filtered result of strong beam
%distance = distance-distance(1);% Set the reference distance to the start of the data
Distance_sig = TBLHC_dbscan(:,1);
height_sig = TBLHC_dbscan(:,2);

%% Display the DBSCAN result of the strong beam
figure()
scatter(TBLHC1(:,1)-TBLHC1(1,1),TBLHC1(:,2),'.');
hold on
scatter(Distance_sig,height_sig,'.');

%% Estimate the along-track slope and noise rate of the strong beam
acc_distance = 20;%accmulated distance
acc_step = 10;%step size
total_window = 1100;%total window width
estimate_max_slope = 50;%estimated maximum slope
threshold_window1 = 50;%the safety range reserved from the terrain profile，+-~m
threshold_window2 = 30;%the safety range reserved from the top and bottom，+-~m
height_window = (total_window -  (threshold_window1+threshold_window2)*2 - tand(estimate_max_slope)*acc_distance)/2;%available window width
segment_pos = (PE_distance(1) + acc_distance/2):acc_step:(PE_distance(end) - acc_distance/2);%split the data into small segments
segment_num = length(segment_pos);
noise_count = zeros(1,segment_num);
slope_seg = zeros(1,segment_num);
slope_conf = zeros(1,segment_num);
index_noise=PE_distance*0;
error_record = segment_pos*0;
tic
for i = 1:segment_num
    %%estimate the noise rate
    index1 = ((PE_distance >= segment_pos(i) - acc_distance/2)...
        .*(PE_distance <= segment_pos(i) + acc_distance/2))>0;
    segment_distance = PE_distance(index1);
    segment_height = PE_height(index1);
    bottom = min(segment_height);    
    top = max(segment_height);
    
    noise_count1 = sum(((segment_height >= bottom + threshold_window2).*...
        (segment_height <= bottom + threshold_window2 + height_window))>0);
    noise_count2 = sum(((segment_height <= top - threshold_window2).*...
        (segment_height >= top - threshold_window2 - height_window))>0);
    noise_count(i) = noise_count1 + noise_count2;
    
    bottom_down = bottom + threshold_window2;
    bottom_up =bottom + threshold_window2 + height_window; 
    top_down = top - threshold_window2 - height_window;
    top_up = top - threshold_window2;
    
    index_noise = index_noise + (((PE_height >= bottom_down).*...
        (PE_height <= bottom_up)) +...
        ((PE_height <= top_up).*...
        (PE_height >=  top_down))).*index1;
    
    
    %%estimate the along-track slope
    index2 = ((Distance_sig >= segment_pos(i) - acc_distance/2)...
        .*(Distance_sig <= segment_pos(i) + acc_distance/2))>0;
    segment_sig_distance = double(Distance_sig(index2));
    segment_sig_height = double(height_sig(index2));
    if length(segment_sig_distance) <=2
        error_record(i) = 1;%recording failed fitting
    end
    [p,s]=polyfit(segment_sig_distance,segment_sig_height,1);%linear fitting
    slope_conf(i) = s.normr/sqrt(length(segment_sig_distance));
    slope_seg(i) = atand(p(1));   
end
toc
noise_rate = noise_count./(height_window*2/0.15*1E-9)/(1E6)/(acc_distance/7.6123E3/1E-4);
% error_record = error_record>0;
% slope_seg(error_record) = [];%remove failed fitting result
% segment_pos(error_record) = [];
% slope_conf(error_record) = [];
% noise_rate(error_record) = [];

figure()%Display the noise rate and along-track slope
subplot(2,1,1)
plot(segment_pos/1E3,noise_rate);
legend('Background Noise rate (MHz)');
grid on
subplot(2,1,2)
bar(segment_pos/1E3,slope_seg);
legend('Along-track Slope (°)');
grid on

TBLHC_step_filterted = [segment_pos',noise_rate',slope_seg',slope_conf'];
%% Fitting the slope-noise relationship
%%ubac slope
slope_step =5;%step size of slope
slope_seg_ubac = 0:slope_step:50;
mean_noise_seg_ubac = zeros(length(slope_seg_ubac),1);
std_noise_seg_ubac = zeros(length(slope_seg_ubac),1);
for i = 1:length(slope_seg_ubac)
    start = slope_seg_ubac(i) - slope_step/2;
    stop = slope_seg_ubac(i) + slope_step/2;
    index = (TBLHC_step_filterted(:,3) >= 0).*(TBLHC_step_filterted(:,3) >= start).*(TBLHC_step_filterted(:,3) < stop)>0;
    mean_noise_seg_ubac(i) = mean(TBLHC_step_filterted(index,2));
    std_noise_seg_ubac(i) = std(TBLHC_step_filterted(index,2));
end


%%adret slope
slope_step =5;%step size of slope
slope_seg_adret = 0:-slope_step:-50;
mean_noise_seg_adret = zeros(length(slope_seg_adret),1);
std_noise_seg_adret = zeros(length(slope_seg_adret),1);
for i = 1:length(slope_seg_adret)
    start = slope_seg_adret(i) - slope_step/2;
    stop = slope_seg_adret(i) + slope_step/2;
    index = (TBLHC_step_filterted(:,3) < 0).*(TBLHC_step_filterted(:,3) >= start).*(TBLHC_step_filterted(:,3) < stop)>0;
    mean_noise_seg_adret(i) = mean(TBLHC_step_filterted(index,2));
    std_noise_seg_adret(i) = std(TBLHC_step_filterted(index,2));
end

[fitresult_adret, gof_adret] = createFit(slope_seg_adret, mean_noise_seg_adret);
adret_a = fitresult_adret.a;
adret_b = fitresult_adret.b;
adret_c = fitresult_adret.c;
adret_d = fitresult_adret.d;
cal_noise_adret = adret_a*cosd(adret_b*slope_seg_adret+adret_c) + adret_d;

[fitresult_ubac, gof_ubac] = createFit(slope_seg_ubac, mean_noise_seg_ubac);
ubac_a = fitresult_ubac.a;
ubac_b = fitresult_ubac.b;
ubac_c = fitresult_ubac.c;
ubac_d = fitresult_ubac.d;
cal_noise_ubac = ubac_a*cosd(ubac_b*slope_seg_ubac+ubac_c) + ubac_d;

%% fitting noise-slope relationship
%%adret
noise_step =0.2;%noise step
noise_seg_adret = 1:noise_step:4;
mean_slope_seg_adret = zeros(length(noise_seg_adret),1);
std_slope_seg_adret = zeros(length(noise_seg_adret),1);
for i = 1:length(noise_seg_adret)
    start = noise_seg_adret(i) - noise_step/2;
    stop = noise_seg_adret(i) + noise_step/2;
    index = (TBLHC_step_filterted(:,2) >= start).*(TBLHC_step_filterted(:,2) < stop).*(TBLHC_step_filterted(:,3) < 0)>0;
    mean_slope_seg_adret(i) = mean(TBLHC_step_filterted(index,3));
    std_slope_seg_adret(i) = std(TBLHC_step_filterted(index,3));
end
%%ubac
noise_step =0.2;%noise step
noise_seg_ubac = 1.6:noise_step:4;
mean_slope_seg_ubac = zeros(length(noise_seg_ubac),1);
std_slope_seg_ubac = zeros(length(noise_seg_ubac),1);
for i = 1:length(noise_seg_ubac)
    start = noise_seg_ubac(i) - noise_step/2;
    stop = noise_seg_ubac(i) + noise_step/2;
    index = (TBLHC_step_filterted(:,2) >= start).*(TBLHC_step_filterted(:,2) < stop).*(TBLHC_step_filterted(:,3) >= 0)>0;
    mean_slope_seg_ubac(i) = mean(TBLHC_step_filterted(index,3));
    std_slope_seg_ubac(i) = std(TBLHC_step_filterted(index,3));
end

%% cosine function fitting
% [fitresult_ns_adret, gof] = createFit2(mean_slope_seg_adret,noise_seg_adret);
% ns_adret_a = fitresult_ns_adret.a;
% ns_adret_b = fitresult_ns_adret.b;
% ns_adret_c = fitresult_ns_adret.c;
% ns_adret_d = fitresult_ns_adret.d;
% cal_slope_adret = (-acosd((noise_seg_adret - ns_adret_d)/ns_adret_a) - ns_adret_c)/ns_adret_b;
% [fitresult_ns_ubac, gof] = createFit2(mean_slope_seg_ubac, noise_seg_ubac);
% ns_ubac_a = fitresult_ns_ubac.a;
% ns_ubac_b = fitresult_ns_ubac.b;
% ns_ubac_c = fitresult_ns_ubac.c;
% ns_ubac_d = fitresult_ns_ubac.d;
% cal_slope_ubac = (-acosd((noise_seg_ubac - ns_ubac_d)/ns_ubac_a) - ns_ubac_c)/ns_ubac_b;

%% Polynomial fitting
%%迎坡
[fitresult_ns_adret, gof] = createFit_polynomial(noise_seg_adret, mean_slope_seg_adret);
ns_adret_p1 = fitresult_ns_adret.p1;
ns_adret_p2 = fitresult_ns_adret.p2;
ns_adret_p3 = fitresult_ns_adret.p3;
ns_adret_p4 = fitresult_ns_adret.p4;
cal_slope_adret = ns_adret_p1*noise_seg_adret.^3+ns_adret_p2*noise_seg_adret.^2 ...
+ns_adret_p3*noise_seg_adret + ns_adret_p4;

%%背坡
[fitresult_ns_ubac, gof] = createFit_polynomial(noise_seg_ubac, mean_slope_seg_ubac);
ns_ubac_p1 = fitresult_ns_ubac.p1;
ns_ubac_p2 = fitresult_ns_ubac.p2;
ns_ubac_p3 = fitresult_ns_ubac.p3;
ns_ubac_p4 = fitresult_ns_ubac.p4;
cal_slope_ubac = ns_ubac_p1*noise_seg_ubac.^3+ns_ubac_p2*noise_seg_ubac.^2 ...
+ns_ubac_p3*noise_seg_ubac + ns_ubac_p4;

%% load weak beam data
% load('TBLHC-2l.mat');
TBLHC2 = TBLHC_2th_weak;
TBLHC2(:,1) = TBLHC2(:,1) - TBLHC2(1,1);  

%% split data into small segment
point_num=25000;% filter point number each time
al_num=ceil(size(TBLHC2,1)/point_num);%valid num points
for kk=1:al_num
    if kk==al_num
%         break;
        d{al_num}=TBLHC2((kk-1)*point_num+1:end,1);
        h{al_num}=TBLHC2((kk-1)*point_num+1:end,2);
        b{al_num}=TBLHC2((kk-1)*point_num+1:end,3);
        l{al_num}=TBLHC2((kk-1)*point_num+1:end,4);
    else
        d{kk}=TBLHC2((kk-1)*point_num+1:kk*point_num,1);
        h{kk}=TBLHC2((kk-1)*point_num+1:kk*point_num,2);
        b{kk}=TBLHC2((kk-1)*point_num+1:kk*point_num,3);%纬度
        l{kk}=TBLHC2((kk-1)*point_num+1:kk*point_num,4);%经度
    end
end
%% mDBSCAN
TBLHC_dbscannew=[];
TBLHC_dbscanold=[];
TBLHC_filtered = cell(al_num,1);
classify_rec = cell(al_num,1);
tic
parfor kk=1:al_num
    [classify_all,classify,dist_r,elev_r,lat_r,lon_r,slope_r,noise_r,minpts_r,snr_r]=filter_DBSCANslope(d{kk},h{kk},b{kk},l{kk},fitresult_ns_adret,fitresult_ns_ubac,segment_pos,slope_seg);
    TBLHC_dbscannew=[];
    TBLHC_dbscannew=[TBLHC_dbscannew;dist_r,elev_r,lat_r,lon_r,dist_r/7.6123e3,slope_r,noise_r,minpts_r,snr_r,classify];% dist H B L GPSTIME
    TBLHC_filtered{kk} = TBLHC_dbscannew;
    classify_rec{kk} = classify_all;
end
toc

PE_type = 0;
temp2 = zeros(1,10);
for i = 1:al_num
    temp2 = [temp2;TBLHC_filtered{i}];
    PE_type = [PE_type;classify_rec{i}];
end
PE_type(1) = [];
TBLHC2 = [TBLHC2,PE_type];
temp2(1,:) = [];

%% Data filter
temp_filter_input2 = CIFilter(TBLHC2);

%% Filtering result
figure()%滤波效果对比
index_final = find(abs(temp_filter_input2(:,end))>0);
scatter(TBLHC2(:,1),TBLHC2(:,2),'r.');
hold on
scatter(temp_filter_input2(index_final,1),temp_filter_input2(index_final,2),'g.');
hold off;
grid on
ylabel('Elevation (km)')
xlabel('Along-track distance (km)')


