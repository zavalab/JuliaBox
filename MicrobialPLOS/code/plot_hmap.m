close all
clear
clc

set(0,'defaultAxesFontSize',14)
speciesOrder={'BH','CA','BU','PC','BO','BV','BT','EL','FP','CH','DP','ER'};
figurePath='../figure/';

fg=figure(1);fg.Position(3:4)=[600 400];
param=reshape(csvread('../output/standard/param.csv'),13,12)';
cmap =[ ones(32,1) linspace(0,1,32)' linspace(0,1,32)';
        linspace(1,0,32)' linspace(1,0,32)' ones(32,1)];
h=heatmap(param');
h.YData = {'-',speciesOrder{:}};
h.XData = speciesOrder;
h.ColorLimits = [-max(abs(param(:))) max(abs(param(:)))];
h.Colormap=cmap;
h.FontSize=20;
h.CellLabelColor='none';
saveas(fg,[figurePath 'hmap_param.eps'],'epsc')

fg=figure(1);fg.Position(3:4)=[600 400];
param=reshape(csvread('../output/mle/param.csv'),13,12)';
cmap =[ ones(32,1) linspace(0,1,32)' linspace(0,1,32)';
        linspace(1,0,32)' linspace(1,0,32)' ones(32,1)];
h=heatmap(param');
h.YData = {'-',speciesOrder{:}};
h.XData = speciesOrder;
h.ColorLimits = [-max(abs(param(:))) max(abs(param(:)))];
h.Colormap=cmap;
h.FontSize=20;
h.CellLabelColor='none';
saveas(fg,[figurePath 'hmap_param_none.eps'],'epsc')

fg=figure(1);fg.Position(3:4)=[600 400];
param=reshape(csvread('../output/L1prior/param.csv'),13,12)';
cmap =[ ones(32,1) linspace(0,1,32)' linspace(0,1,32)';
        linspace(1,0,32)' linspace(1,0,32)' ones(32,1)];
h=heatmap(param');
h.YData = {'-',speciesOrder{:}};
h.XData = speciesOrder;
h.ColorLimits = [-max(abs(param(:))) max(abs(param(:)))];
h.Colormap=cmap;
h.FontSize=20;
h.CellLabelColor='none';
saveas(fg,[figurePath 'hmap_param_L1.eps'],'epsc')

