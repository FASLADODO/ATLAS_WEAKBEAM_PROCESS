function [ temp_filter_input2 ] = SSDF( temp_filter_input2,Neighbours,slope_e,slope)
%Small-Scale Spatial Density Filter
%   此处显示详细说明
num = sum(Neighbours);
seq = find(Neighbours>0);
a = 17/2;
% b = a*abs(tand(slope_e));
b=a/2;
R = [cosd(slope),-sind(slope); sind(slope),cosd(slope)];
noise_aver = mean(temp_filter_input2(:,end-1));
for i = 1:num
    x0=temp_filter_input2(seq(i),1);
    y0=temp_filter_input2(seq(i),2);
    TH = temp_filter_input2(seq(i),end-1)*(2*b/0.15)*1E-3*(17/0.7)*pi/4;%threshold
    if TH==0
        TH = noise_aver*(2*b/0.15)*1E-3*(17/0.7)*pi/4;
    end
    c_orig = [temp_filter_input2(:,1)-x0,temp_filter_input2(:,2)-y0];
    c_rot = c_orig*R;
    density = sum((c_rot(:,1).^2/a^2 + c_rot(:,2).^2/b^2)<=1);%spatial density
    if density >= 3*TH
        temp_filter_input2(seq(i),end) = -2;
    end
end

end

