clc;
clear all;
close all;

% Reading Image
img = imread('TP_better.png'); axis off;
img_gray = rgb2gray(img);
% Reading Node Locations

budget = '1.82e8';
epsilon = '1';
str_0 = sprintf('technologies_sited_%s.csv', epsilon);
%str_0 = sprintf('node_matrix.csv');

fidi = fopen(str_0);
d = textscan(fidi, '%s %f %f %s %f', 'Delimiter',',', 'HeaderLines',1);
fclose(fidi);

node = d{1};
lat = cell2mat(d(2));
long = cell2mat(d(3));
tech_type = d{4};
tech_size = cell2mat(d(5));

alpha = 1.5757022453;
hold on;
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

lat_plot = a + b*lat_rot;
long_plot = c + d*long_rot;

fig = figure(1);
%imshow(img_gray,'InitialMagnification','fit'); %image(image)
imshow(img_gray);
hold on;

location_plot = [round(lat_plot), round(long_plot)];

r = ones(length(lat),1);
g = ones(length(lat),1);
b = ones(length(lat),1);


for i = 1:length(lat)
   r(i) = img(round(long_plot(i)),round(lat_plot(i)),1);
   g(i) = img(round(long_plot(i)),round(lat_plot(i)), 2);
   b(i) = img(round(long_plot(i)),round(lat_plot(i)), 3);
end

pixel_value = [r g b];

img_scaled = double(img/255);
hsl = rgb2hsl(img_scaled);

h = ones(length(lat),1);
s = ones(length(lat),1);
l = ones(length(lat),1);


for i = 1:length(lat)
   h(i) = hsl(round(long_plot(i)),round(lat_plot(i)), 1)*360;
   s(i) = hsl(round(long_plot(i)),round(lat_plot(i)), 2);
   l(i) = hsl(round(long_plot(i)),round(lat_plot(i)), 3)*100;
end

conc = ones(length(lat),1);

for j = 1:length(lat)
   if 230 < h(j) && h(j) <= 250  
       conc(j) = 1;
       
   elseif 210 < h(j) && h(j) <= 230
       conc(j) = 2;

   elseif 180 < h(j) && h(j) <= 210
       conc(j) = 3;
   elseif 150 < h(j) && h(j) <= 180
       conc(j) = 4;
   elseif 120 < h(j) && h(j) <= 150
       conc(j) = 5;       
   elseif 100 < h(j) && h(j) <= 120
       conc(j) = 6;    
   elseif 80 < h(j) && h(j) <= 100
       conc(j) = 7;    
   elseif 70 < h(j) && h(j) <= 80
       conc(j) = 8;    
   elseif 60 < h(j) && h(j) <= 70
       conc(j) = 9;    
   elseif 0 < h(j) && h(j) <= 60
       conc(j) = 10;    
   end
     
end

% scatter(lat_plot, long_plot)
%scatter_patches(lat_plot, long_plot, (h+1)/10, 'r','FaceAlpha',0.6,'EdgeColor','none');

%fig = figure;
%hs = scatter_patches(lat_plot, long_plot, 20,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','none','MarkerSize',1);
%scatter(lat_plot, long_plot,100,'k', 'filled')
%scatter(lat_plot, long_plot,100*tech_size/max(tech_size),'k', 'filled')
hold on;
%gscatter(lat_plot, long_plot, tech_type,'bryp',[],20*tech_size/max(tech_size));
axis off;
%set(hs,'MarkerFaceColor','r');
%alpha(hs,.5);
%print(fig,'20mil_1','-dpng')

%% Plotting flows
% f(1) = fopen('flow_p1_results.csv');
% f(2) = fopen('flow_p2_results.csv');
% f(3) = fopen('flow_p3_results.csv');

%f = [f1, f2, f3];
plist = [{'p1'}, {'p2'},{'p3'},{'p4'},{'p5'},{'p6'},{'p7'},{'p8'},{'p9'},{'p10'},{'p11'},{'p12'},{'p13'},{'p14'},{'p15'},{'p16'}];    %list of products

nCols = 444;                        %no. of nodes
format = ['%s' repmat(' %f', [1 nCols])];   %format of data in the result file, remat helps in specifying data for 100 columns

%link(1,1,1:3) = 0;

prod = 10;
flow_tot = zeros(prod,1);
dist_tot = zeros(prod,1);
for v=1:prod
    
    str = sprintf('flow_%s_results_%s.csv',char(plist(v)), epsilon);
    f = fopen(str);
    d = textscan(f, format, 'Delimiter',',', 'HeaderLines',1);
    fclose(f);

%link = [];
    q = 1;  
    for k= 2:(nCols+1)
        stor = cell2mat(d(k)); 
   
        for i = 1:nCols
              if stor(i) >= 1
                  link(q,1,v) = i;
                  link(q,2,v) = k-1;
                  q = q+1;
              end
            flow_tot(v) = flow_tot(v) + stor(i);
        end   
        
    end
end
% fclose(f1);
% fclose(f2);
% fclose(f3);
%%%%
 ref = fopen('node_matrix.csv');
 r = textscan(ref, '%s %s %f %f %f', 'Delimiter',',', 'HeaderLines',1);
 fclose(ref);

lat_r = cell2mat(r(3));
long_r = cell2mat(r(4));

j=1;
for i = 1:prod
    
    try 
    for k =1:length(link(:,1,i))
        if link(k,1,i)~=0 && link(k,2,i) ~=0
            lat_p(j,1,i) = lat_r(link(k,1,i));
            long_p(j,1,i) = long_r(link(k,1,i));
            %j = j+1;
            lat_p(j,2,i) = lat_r(link(k,2,i));
            long_p(j,2,i) = long_r(link(k,2,i));
            j = j+1;
        end
    end
    end
    
    %trial(~any(trial,2),:) = [];  %% to delete rows with 0 entries
end

try
lat_p1 = lat_p(:,:,1);
long_p1 = long_p(:,:,1);
end

try 
lat_p2 = lat_p(:,:,2);
long_p2 = long_p(:,:,2);
end

try
lat_p3 = lat_p(:,:,3);
long_p3 = long_p(:,:,3);
end

try
lat_p4 = lat_p(:,:,4);
long_p4 = long_p(:,:,4);
end

try
lat_p5 = lat_p(:,:,5);
long_p5 = long_p(:,:,5);
end

try
lat_p6 = lat_p(:,:,6);
long_p6 = long_p(:,:,6);
end

try
lat_p7 = lat_p(:,:,7);
long_p7 = long_p(:,:,7);
end

try
lat_p8 = lat_p(:,:,8);
long_p8 = long_p(:,:,8);
end

try
lat_p9 = lat_p(:,:,9);
long_p9 = long_p(:,:,9);
end

try
lat_p10 = lat_p(:,:,10);
long_p10 = long_p(:,:,10);
end

try
lat_p11 = lat_p(:,:,11);
long_p11 = long_p(:,:,11);
end

try
lat_p12 = lat_p(:,:,12);
long_p12 = long_p(:,:,12);
end

try
lat_p13 = lat_p(:,:,13);
long_p13 = long_p(:,:,13);
end

try
lat_p14 = lat_p(:,:,14);
long_p14 = long_p(:,:,14);
end

try
lat_p15 = lat_p(:,:,15);
long_p15 = long_p(:,:,15);
end

try
lat_p16 = lat_p(:,:,16);
long_p16 = long_p(:,:,16);
end
%% Removing 0 elements

try
lat_p2(~any(lat_p2,2),:) = [];
long_p2(~any(long_p2,2),:) = [];

[pixel_2x, pixel_2y] = tp_map_better(lat_p2, long_p2); % using the function defined for calculation of pixel co-ordinates
h2 = plot(pixel_2x', pixel_2y', 'b', 'LineWidth', 7.5);

catch
    fprintf('No links for prodcut p2');
end
%p2:sludge
try 
lat_p3(~any(lat_p3,2),:) = [];
long_p3(~any(long_p3,2),:) = [];

[pixel_3x, pixel_3y] = tp_map_better(lat_p3, long_p3); % using the function defined for calculation of pixel co-ordinates
h3 = plot(pixel_3x', pixel_3y', 'r', 'LineWidth', 3);

catch
    fprintf('No links for prodcut p3 \n');
end
%p3:foodwaste
try
lat_p1(~any(lat_p1,2),:) = [];
long_p1(~any(long_p1,2),:) = [];

[pixel_1x, pixel_1y] = tp_map_better(lat_p1, long_p1); % using the function defined for calculation of pixel co-ordinates
h1 = plot(pixel_1x', pixel_1y',  'g', 'LineWidth', 7.5);
hold on;

catch
    fprintf('No links for prodcut p1');

end


try
lat_p4(~any(lat_p4,2),:) = [];
long_p4(~any(long_p4,2),:) = [];

[pixel_4x, pixel_4y] = tp_map_better(lat_p4, long_p4); % using the function defined for calculation of pixel co-ordinates
h4 = plot(pixel_4x', pixel_4y', 'y', 'LineWidth', 3);

catch
    fprintf('No links for prodcut p4');
end
%p5:c8

try
lat_p5(~any(lat_p5,2),:) = [];
long_p5(~any(long_p5,2),:) = [];

[pixel_5x, pixel_5y] = tp_map_better(lat_p5, long_p5); % using the function defined for calculation of pixel co-ordinates
h5 = plot(pixel_5x', pixel_5y', 'y', 'LineWidth', 3);

catch
    fprintf('No links for prodcut p5');
end
%p6:digestate

try
lat_p6(~any(lat_p6,2),:) = [];
long_p6(~any(long_p6,2),:) = [];

[pixel_6x, pixel_6y] = tp_map_better(lat_p6, long_p6); % using the function defined for calculation of pixel co-ordinates
h6 = plot(pixel_6x', pixel_6y', 'k', 'LineWidth', 3);

catch
    fprintf('No links for prodcut p6');
end

try
lat_p7(~any(lat_p7,2),:) = [];
long_p7(~any(long_p7,2),:) = [];

[pixel_7x, pixel_7y] = tp_map_better(lat_p7, long_p7); % using the function defined for calculation of pixel co-ordinates
h7 = plot(pixel_7x', pixel_7y', 'k', 'LineWidth', 3);

catch
    fprintf('No links for prodcut p7');
end

try
lat_p8(~any(lat_p8,2),:) = [];
long_p8(~any(long_p8,2),:) = [];

[pixel_8x, pixel_8y] = tp_map_better(lat_p8, long_p8); % using the function defined for calculation of pixel co-ordinates
h8 = plot(pixel_8x', pixel_8y', 'k', 'LineWidth', 3);

catch
    fprintf('No links for prodcut p6');
end

try
lat_p9(~any(lat_p9,2),:) = [];
long_p9(~any(long_p9,2),:) = [];

[pixel_9x, pixel_9y] = tp_map_better(lat_p9, long_p9); % using the function defined for calculation of pixel co-ordinates
h9 = plot(pixel_9x', pixel_9y', 'k', 'LineWidth', 3);

catch
    fprintf('No links for prodcut p9');
end

try
lat_p10(~any(lat_p10,2),:) = [];
long_p10(~any(long_p10,2),:) = [];

[pixel_10x, pixel_10y] = tp_map_better(lat_p10, long_p10); % using the function defined for calculation of pixel co-ordinates
h10 = plot(pixel_10x', pixel_10y', 'k', 'LineWidth', 3);

catch
    fprintf('No links for prodcut p10');
end

try
lat_p11(~any(lat_p11,2),:) = [];
long_p11(~any(long_p11,2),:) = [];

[pixel_11x, pixel_11y] = tp_map_better(lat_p11, long_p11); % using the function defined for calculation of pixel co-ordinates
h10 = plot(pixel_11x', pixel_11y', 'k', 'LineWidth', 3);

catch
    fprintf('No links for prodcut p11');
end

try
lat_p12(~any(lat_p12,2),:) = [];
long_p12(~any(long_p12,2),:) = [];

[pixel_12x, pixel_12y] = tp_map_better(lat_p12, long_p12); % using the function defined for calculation of pixel co-ordinates
h12 = plot(pixel_12x', pixel_12y', 'k', 'LineWidth', 3);

catch
    fprintf('No links for prodcut p12');
end

try
lat_p13(~any(lat_p13,2),:) = [];
long_p13(~any(long_p13,2),:) = [];

[pixel_13x, pixel_13y] = tp_map_better(lat_p13, long_p13); % using the function defined for calculation of pixel co-ordinates
h13 = plot(pixel_13x', pixel_13y', 'k', 'LineWidth', 3);

catch
    fprintf('No links for prodcut p13');
end

try
lat_p14(~any(lat_p14,2),:) = [];
long_p14(~any(long_p14,2),:) = [];

[pixel_14x, pixel_14y] = tp_map_better(lat_p14, long_p14); % using the function defined for calculation of pixel co-ordinates
h14 = plot(pixel_14x', pixel_14y', 'k', 'LineWidth', 3);

catch
    fprintf('No links for prodcut p14');
end

try
lat_p15(~any(lat_p15,2),:) = [];
long_p15(~any(long_p15,2),:) = [];

[pixel_15x, pixel_15y] = tp_map_better(lat_p15, long_p15); % using the function defined for calculation of pixel co-ordinates
h15 = plot(pixel_15x', pixel_15y', 'k', 'LineWidth', 3);

catch
    fprintf('No links for prodcut p15');
end

try
lat_p16(~any(lat_p16,2),:) = [];
long_p16(~any(long_p16,2),:) = [];

[pixel_16x, pixel_16y] = tp_map_better(lat_p16, long_p16); % using the function defined for calculation of pixel co-ordinates
h16 = plot(pixel_16x', pixel_16y', 'k', 'LineWidth', 3);

catch
    fprintf('No links for prodcut p16');
end

[pixel_rx, pixel_ry] = tp_map_better(lat_r, long_r);


%scatter(lat_plot, long_plot,100,'k', 'filled')

index_1 = find(strcmp(tech_type, 'tA1')+strcmp(tech_type, 'tA2')+strcmp(tech_type, 'tA3')+strcmp(tech_type, 'tA4'));
index_2 = find(strcmp(tech_type, 'tB1')+strcmp(tech_type, 'tB2')+strcmp(tech_type, 'tB3')+strcmp(tech_type, 'tB4'));
index_3 = find(strcmp(tech_type, 'tC1')+strcmp(tech_type, 'tC2')+strcmp(tech_type, 'tC3')+strcmp(tech_type, 'tC4'));
index_4 = find(strcmp(tech_type, 'tD1')+strcmp(tech_type, 'tD2')+strcmp(tech_type, 'tD3')+strcmp(tech_type, 'tD4'));
index_5 = find(strcmp(tech_type, 'tE1')+strcmp(tech_type, 'tE2')+strcmp(tech_type, 'tE3')+strcmp(tech_type, 'tE4'));
index_6 = find(strcmp(tech_type, 'tF1')+strcmp(tech_type, 'tF2')+strcmp(tech_type, 'tF3')+strcmp(tech_type, 'tF4')+strcmp(tech_type,'tG1')+strcmp(tech_type, 'tG2')+strcmp(tech_type, 'tG3')+strcmp(tech_type, 'tG4')+strcmp(tech_type,'tH1')+strcmp(tech_type, 'tH2')+strcmp(tech_type, 'tH3')+strcmp(tech_type, 'tH4')+strcmp(tech_type,'tI1')+strcmp(tech_type, 'tI2')+strcmp(tech_type, 'tI3')+strcmp(tech_type, 'tI4')+strcmp(tech_type,'tJ1')+strcmp(tech_type, 'tJ2')+strcmp(tech_type, 'tJ3')+strcmp(tech_type, 'tJ4')+strcmp(tech_type,'tK1')+strcmp(tech_type, 'tK2')+strcmp(tech_type, 'tK3')+strcmp(tech_type, 'tK4')+strcmp(tech_type,'tL1')+strcmp(tech_type, 'tL2')+strcmp(tech_type, 'tL3')+strcmp(tech_type, 'tL4')+strcmp(tech_type,'tM1')+strcmp(tech_type, 'tM2')+strcmp(tech_type, 'tM3')+strcmp(tech_type, 'tM4')+strcmp(tech_type,'tN1')+strcmp(tech_type, 'tN2')+strcmp(tech_type, 'tN3')+strcmp(tech_type, 'tN4')+strcmp(tech_type,'tO1')+strcmp(tech_type, 'tO2')+strcmp(tech_type, 'tO3')+strcmp(tech_type, 'tO4'));

for i = 1:length(index_6)
    scatter(lat_plot(index_6(i)), long_plot(index_6(i)),30,[1 0 0], 'filled')
end 

for i = 1:length(index_1)
    scatter(lat_plot(index_1(i)), long_plot(index_1(i)),15,'g', 'filled')
end 

for i = 1:length(index_2)
    scatter(lat_plot(index_2(i)), long_plot(index_2(i)),15,'w', 'filled')
end 

for i = 1:length(index_3)
    scatter(lat_plot(index_3(i)), long_plot(index_3(i)),15,[1 1 0], 'filled')
end 

for i = 1:length(index_4)
    scatter(lat_plot(index_4(i)), long_plot(index_4(i)),15,'g', 'filled')
end 

for i = 1:length(index_5)
    scatter(lat_plot(index_5(i)), long_plot(index_5(i)),15,'w', 'filled')
end 


h11 = scatter(pixel_rx(1:72), pixel_ry(1:72),6, 'k', 'filled');
h12 = scatter(pixel_rx(73:444), pixel_ry(73:444),6, 'b', 'filled');

%scatter(pixel_rx(nCols), pixel_ry(nCols),100, 'g', 'filled')
%scatter(pixel_rx(nCols), pixel_ry(nCols), 'g', 'filled','MarkerEdgeColor','y')
scatter(pixel_rx(nCols), pixel_ry(nCols), 30, [0.9961 0.9687 0.0742], 'filled')
%%%%

str_p = sprintf('%sk_%s', budget, epsilon);


%legend([h2 h4], {'Struvite','farms'});

print(fig, str_p,'-depsc')   % prints a color version in eps format

% finding indexes

% Output: 70cm * 88.47cm

