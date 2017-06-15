clc;
clear all;
close all;

% Reading Image
img = imread('TP_better.png'); axis off;
img_gray = rgb2gray(img);
% Reading Node Locations

budget = '5k';
epsilon = '0.5';
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
plist = [{'p1'}, {'p2'},{'p3'}];    %list of products

nCols = 101;                        %no. of nodes
format = ['%s' repmat(' %f', [1 nCols])];   %format of data in the result file, remat helps in specifying data for 100 columns

%link(1,1,1:3) = 0;

prod = 3;
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

%% Removing 0 elements

try
lat_p2(~any(lat_p2,2),:) = [];
long_p2(~any(long_p2,2),:) = [];

[pixel_2x, pixel_2y] = tp_map_better(lat_p2, long_p2); % using the function defined for calculation of pixel co-ordinates
h2 = plot(pixel_2x', pixel_2y', 'y', 'LineWidth', 1.5);

catch
    fprintf('No links for prodcut p2');
end

try 
lat_p3(~any(lat_p3,2),:) = [];
long_p3(~any(long_p3,2),:) = [];

[pixel_3x, pixel_3y] = tp_map_better(lat_p3, long_p3); % using the function defined for calculation of pixel co-ordinates
h3 = plot(pixel_3x', pixel_3y', 'k', 'LineWidth', 1.5);
%legend([h3],{'Digestate'})

catch
    fprintf('No links for prodcut p3 \n');
end

try
lat_p1(~any(lat_p1,2),:) = [];
long_p1(~any(long_p1,2),:) = [];

[pixel_1x, pixel_1y] = tp_map_better(lat_p1, long_p1); % using the function defined for calculation of pixel co-ordinates
h1 = plot(pixel_1x', pixel_1y',  'b', 'LineWidth', 1.5);
hold on;

catch
    fprintf('No links for prodcut p1');

end







[pixel_rx, pixel_ry] = tp_map_better(lat_r, long_r);


%scatter(lat_plot, long_plot,100,'k', 'filled')

index_1 = find(strcmp(tech_type, 't1'));
index_2 = find(strcmp(tech_type, 't2'));
index_3 = find(strcmp(tech_type, 't3'));
index_4 = find(strcmp(tech_type, 't4'));
index_5 = find(strcmp(tech_type, 't5'));
index_6 = find(strcmp(tech_type, 't6'));
index_7 = find(strcmp(tech_type, 't7'));
index_8 = find(strcmp(tech_type, 't8'));
index_9 = find(strcmp(tech_type, 't9'));

for i = 1:length(index_1)
    scatter(lat_plot(index_1(i)), long_plot(index_1(i)),100,'y', 'filled')
end 

for i = 1:length(index_2)
    scatter(lat_plot(index_2(i)), long_plot(index_2(i)),100,'y', 'filled')
end 

for i = 1:length(index_3)
    scatter(lat_plot(index_3(i)), long_plot(index_3(i)),100,'y', 'filled')
end 

for i = 1:length(index_4)
    scatter(lat_plot(index_4(i)), long_plot(index_4(i)),100,[0 0.5 0], 'filled')
end 

for i = 1:length(index_5)
    scatter(lat_plot(index_5(i)), long_plot(index_5(i)),100,[0 0.5 0], 'filled')
end 

for i = 1:length(index_6)
    scatter(lat_plot(index_6(i)), long_plot(index_6(i)),100,[0 0.5 0], 'filled')
end 

for i = 1:length(index_7)
    scatter(lat_plot(index_7(i)), long_plot(index_7(i)),100,'k', 'filled')
end 

for i = 1:length(index_8)
    scatter(lat_plot(index_8(i)), long_plot(index_8(i)),100,'k', 'filled')
end 

for i = 1:length(index_9)
    scatter(lat_plot(index_9(i)), long_plot(index_9(i)),100,'k', 'filled')
end 

scatter(pixel_rx, pixel_ry,'r', 'filled');

%scatter(pixel_rx(nCols), pixel_ry(nCols),100, 'g', 'filled')
scatter(pixel_rx(nCols), pixel_ry(nCols), 'g', 'filled','MarkerEdgeColor','y')

%%%%

str_p = sprintf('alpha_%s', epsilon);


%legend([h1 h2 h3], {'Manure','Struvite','Digestate'});

print(fig, str_p,'-dpng')

% finding indexes


