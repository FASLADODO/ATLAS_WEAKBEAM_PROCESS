function [ r_e] = DrawEllipse( x1,y1,x2,y2,slope,A,B)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
pos_x = (x1+x2)/2;
pos_y = (y1+y2)/2;
% A = sqrt((x2-x1)^2+(y2-y1)^2)/2;
% B = A/2;
theta_grid = linspace(0,2*pi);
e_x = A*cos(theta_grid);
e_y = B*sin(theta_grid);
phi = slope;
R=  [ cosd(phi),sind(phi); -sind(phi),cosd(phi)];
r_e = [e_x;e_y]'*R;
r_e(:,1) = r_e(:,1) + pos_x;
r_e(:,2) = r_e(:,2) + pos_y;
end

