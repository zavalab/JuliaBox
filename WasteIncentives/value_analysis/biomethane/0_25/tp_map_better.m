function [ pixel_x, pixel_y ] = tp_map( lat, long )
%Plot lat-long data on TP concentration image file
%   Detailed explanation goes here
alpha = 1.5757022453;

lat_ref = 42.508307;	
long_ref = -90.640803;

location = [lat, long];
% 
% latitude = location(:,1);
% longitude = location(:,2);

lat_scaled = lat - lat_ref;
long_scaled = long - long_ref;

lat_rot = lat_scaled*cos(alpha) + long_scaled*sin(alpha);
long_rot = -lat_scaled*sin(alpha) + long_scaled*cos(alpha);

xplot_ref = 420;
yplot_ref = 1263;

a =	436.1742226003;		
b = 192.6414116036;		

c =	1261.9999260322;
d = 272.2654880799;

pixel_x = a + b*lat_rot;
pixel_y = c + d*long_rot;


end

