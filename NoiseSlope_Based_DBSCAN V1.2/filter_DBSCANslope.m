function [classify_all,classify_r,dist_r,elev_r,lat_r,lon_r,slope_r,noise_r,minpts_r,snr_r]=filter_DBSCANslope(dist_initial,elev_initial,lat_initial,lon_initial,fitresult_ns_adret,fitresult_ns_ubac,segment_pos,slope_seg)
Origin_data=[dist_initial elev_initial];
Ns0 = 0.5548;% Coefficient for calculating Average signal photon number
estimate_max_slope = 50;
threshold_window1 = 50;
threshold_window2 = 30;
total_window = 1100;% Must be modified accroding to the range gate of dataset
acc_distance = 20;
Flag=DBSCANslope(Origin_data,fitresult_ns_adret,fitresult_ns_ubac,segment_pos,slope_seg,Ns0,total_window,estimate_max_slope,threshold_window1,threshold_window2,acc_distance);
B=find(Flag(:,5)~=0);
classify_all = Flag(:,5);
%% save data to dist_r elev_result lat_result lon_result
k=1;
for j=1:length(B)
    lat_r(k)=lat_initial(B(j));
    lon_r(k)=lon_initial(B(j));
    elev_r(k)=elev_initial(B(j));
    dist_r(k)=dist_initial(B(j));
    slope_r(k) = Flag(B(j),1);
    noise_r(k) = Flag(B(j),2);
    minpts_r(k) = Flag(B(j),3);
    snr_r(k) = Flag(B(j),4);
    classify_r(k) =Flag(B(j),5);
    k=k+1;
end;
lat_r=lat_r';
lon_r=lon_r';
elev_r=elev_r';
dist_r=dist_r';
slope_r = slope_r';
noise_r = noise_r';
minpts_r = minpts_r';
snr_r = snr_r';
classify_r = classify_r';
