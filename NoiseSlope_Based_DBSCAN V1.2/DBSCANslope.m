%This programm is developed based on the DBSCAN algorithm of Yarpiz team.
% The copyright of the original DBSCAN algorithm belongs to the original
% developers. The copyright statement of the original DBSCAN algorithm is
% as follows:
%
%Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPML110
% Project Title: Implementation of DBSCAN Clustering in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com

%first column Distance, second column elevation
function [IDX,isnoise]=DBSCANslope(X,fitresult_ns_adret,fitresult_ns_ubac,segment_pos,slope_seg,Ns0,total_window,estimate_max_slope,threshold_window1,threshold_window2,acc_distance)
%%
global FLAG
FLAG=0;
n=size(X,1);
IDX=zeros(n,5);%1st column slope,2nd noise,3rd Minpts, 5th classify
visited=false(n,1);
isnoise=false(n,1);

height_window = (total_window -  (threshold_window1+threshold_window2)*2 - tand(estimate_max_slope)*acc_distance)/2;
    for i=1:n
        if ~visited(i)
            visited(i)=true;
            dis_range = 10;
            Dis_range_left = X(i,1)-dis_range;
            Dis_range_right = X(i,1)+dis_range;
            index_dis = X(:,1)>Dis_range_left & X(:,1)<=Dis_range_right;
            X_dis_range = X(index_dis,:);
            Dis_real_left = min(X_dis_range(:,1));
            Dis_real_right = max(X_dis_range(:,1));
            Real_dis_range = Dis_real_right - Dis_real_left;
            altitude_high = max(X_dis_range(:,2));
            altitude_low = min(X_dis_range(:,2));
            noise = sum(((X_dis_range(:,2) >= altitude_low+threshold_window2).*(X_dis_range(:,2) <= altitude_low+threshold_window2+height_window)...
                +(X_dis_range(:,2) <= altitude_high-threshold_window2).*(X_dis_range(:,2) >= altitude_high-threshold_window2-height_window)) > 0);
            norm_noise = noise/((height_window)*2/0.15*1E-9)/(1E6)/(Real_dis_range/7.6123E3/1E-4);
            adret_slope = fitresult_ns_adret.p1*norm_noise.^3 + fitresult_ns_adret.p2*norm_noise.^2 + fitresult_ns_adret.p3*norm_noise + fitresult_ns_adret.p4;
            ubac_slope = fitresult_ns_ubac.p1*norm_noise.^3 + fitresult_ns_ubac.p2*norm_noise.^2 + fitresult_ns_ubac.p3*norm_noise + fitresult_ns_ubac.p4;

            [Neighbors_adret,A1,B1,pulse_sigma1]=RegionQuery(i,adret_slope);
            [Neighbors_ubac,A2,B2,pulse_sigma2]=RegionQuery(i,ubac_slope);
            noise_TH1 = norm_noise*1E6*(B1*2/0.15*1E-9)*(A1*2/7.6123E3/1E-4)*pi/4;
            noise_TH2 = norm_noise*1E6*(B2*2/0.15*1E-9)*(A2*2/7.6123E3/1E-4)*pi/4;
            if numel(Neighbors_adret) >= numel(Neighbors_ubac)
                TH = noise_TH1;
                slope = adret_slope;
                B=B1;
                pulse_sigma = pulse_sigma1;
            else
                TH = noise_TH2;
                slope = ubac_slope;
                B=B2;
                pulse_sigma = pulse_sigma2;
            end
            MinPts = round(3*TH);
            Ns = Ns0*(A1*2/7.6123E3/1E-4)*cosd(slope)*pi/4*4/6;%Estimated average signal photon number
            SNR = 10*log10((Ns0*cosd(slope)+TH*1E6*pulse_sigma*2.355)/(TH*1E6*pulse_sigma*2.355));
            if MinPts < Ns+TH
                MinPts = Ns+TH;
            end
            if MinPts >TH*SNR
                MinPts = TH*SNR;
            end         
             
            if (numel(Neighbors_adret)<MinPts && numel(Neighbors_ubac)<MinPts) 
                isnoise(i)=true;
            end
            if (numel(Neighbors_adret)>=MinPts) || (numel(Neighbors_ubac)>=MinPts)
                if numel(Neighbors_adret) >= numel(Neighbors_ubac)
                    ExpandCluster(i,Neighbors_adret,FLAG);
                else
                    ExpandCluster(i,Neighbors_ubac,FLAG);
                end
            end
        end
    end

    function ExpandCluster(i,Neighbors)
        FLAG = FLAG + 1;
        IDX(i,1)=FLAG;
        k = 1;
        while true
            j = Neighbors(k);
            if ~visited(j)
                visited(j)=true;
                Dis_range_left = X(j,1)-dis_range;
                Dis_range_right = X(j,1)+dis_range;
                index_dis = X(:,1)>Dis_range_left & X(:,1)<Dis_range_right;
                X_dis_range = X(index_dis,:);
                Dis_real_left = min(X_dis_range(:,1));
                Dis_real_right = max(X_dis_range(:,1));
                Real_dis_range = Dis_real_right - Dis_real_left;
                altitude_high = max(X_dis_range(:,2));
                altitude_low = min(X_dis_range(:,2));
                noise = sum(((X_dis_range(:,2) >= altitude_low+threshold_window2).*(X_dis_range(:,2) <= altitude_low+threshold_window2+height_window)...
                    +(X_dis_range(:,2) <= altitude_high-threshold_window2).*(X_dis_range(:,2) >= altitude_high-threshold_window2-height_window)) > 0);
                norm_noise = noise/((height_window)*2/0.15*1E-9)/(1E6)/(Real_dis_range/7.6123E3/1E-4);
                adret_slope = fitresult_ns_adret.p1*norm_noise.^3 + fitresult_ns_adret.p2*norm_noise.^2 + fitresult_ns_adret.p3*norm_noise + fitresult_ns_adret.p4;
                ubac_slope = fitresult_ns_ubac.p1*norm_noise.^3 + fitresult_ns_ubac.p2*norm_noise.^2 + fitresult_ns_ubac.p3*norm_noise + fitresult_ns_ubac.p4;

                [Neighbors_adret,A1,B1]=RegionQuery(j,adret_slope);%return region smaller than radius
                [Neighbors_ubac,A2,B2]=RegionQuery(j,ubac_slope);%return region smaller than radius
                if (numel(Neighbors_adret)<MinPts && numel(Neighbors_ubac)<MinPts) 
                    isnoise(j)=true;
                end
                if (numel(Neighbors_adret)>=MinPts) || (numel(Neighbors_ubac)>=MinPts)
                    if numel(Neighbors_adret)>=numel(Neighbors_ubac)
                        IDX(j,1)=adret_slope;
                        IDX(j,2)=norm_noise;
                        IDX(j,3)=MinPts;
                        IDX(j,4)=SNR;
                        IDX(j,5)=FLAG;
                        Neighbors=[Neighbors; Neighbors_adret];
                    else
                        IDX(j,1)=ubac_slope;
                        IDX(j,2)=norm_noise;
                        IDX(j,3)=MinPts;
                        IDX(j,4)=SNR;
                        IDX(j,5)=FLAG;
                        Neighbors=[Neighbors; Neighbors_ubac];
                    end                    
                end
            end
            if IDX(j,1)==0
                IDX(j,1)=361;%not core
                IDX(j,2)=norm_noise;
                IDX(j,3)=MinPts;
                IDX(j,4)=SNR;
            end

            k = k + 1;
            if k > numel(Neighbors)
                break;
            end
        end
    end

    function [Neighbors,A,B,pulse_sigma]=RegionQuery(i,slope)
        theta = slope;
        Var_land = 0.01;
        pulse_width0 = 1.5E-9;%pulse width,ns
        pulse_sigma0 = pulse_width0/2.355;%RMS PW,ns
        z = 500E3;
        theta_t = 17.5/z/4;       
        c = 2.997E8;
        pulse_sigma = pulse_sigma0^2 + 4*Var_land/c^2 + 4*z^2/c^2.*(tan(theta_t)^4 + tan(theta_t)^2*tand(theta).^2);
        pulse_sigma = sqrt(pulse_sigma);
        B = pulse_sigma*1E9*0.15*4/2;
        A = acc_distance/2;
        TransX = X(:,1) - X(i,1);
        TransY = X(:,2) - X(i,2);
        TransX2 = TransX * cosd(theta) + TransY * sind(theta);
        TransY2 = -TransX * sind(theta) + TransY * cosd(theta);
        dis = power(TransX2,2)/power(A,2) + power(TransY2,2)/power(B,2);
        Neighbors=find(dis<1);
    end
    

end



