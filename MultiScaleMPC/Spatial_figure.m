close all
mkdir figure

clm=colormap;
clm2=colormap('hsv');
colormap default

dFull=csvread('output/Spatial/dFull.csv');
dWaveFull=csvread('output/Spatial/dWaveFull.csv');
dPeakFull=csvread('output/Spatial/dPeakFull.csv');

zFull=csvread('output/Spatial/zFull.csv');
zGS=csvread('output/Spatial/zGS.csv');
zsetFull=csvread('output/Spatial/zsetFull.csv');
zCoarse=csvread('output/Spatial/zCoarse.csv');
zSave1=csvread('output/Spatial/zSave1.csv');
zSave2=csvread('output/Spatial/zSave2.csv');
zSaveSettled=csvread('output/Spatial/zSave10.csv');
errorSave=csvread('output/Spatial/errorSave.csv');

zMax=max([max(zFull(:)),max(zSave1(:)),max(zSave2(:)),max(zCoarse(:))]);
zMin=min([min(zFull(:)),min(zSave1(:)),min(zSave2(:)),min(zCoarse(:))]);
zMaxError=zMax;
zMinError=zMin;
dMax=max(dFull(:));
dMin=min(dFull(:));

figureSize=[18 18];

figure(1)
clf

%%
subplot(3,2,5)
surf(zSaveSettled)
zlabel('State') 
zlim([zMin zMax])
set(gca,'fontsize',16,'Position',[0.08 0.03 0.4 0.30],'TitleFontWeight','normal','TitleFontSizeMultiplier',1.0)
box on

subplot(3,2,6)
surf(zSaveSettled-zFull)
zlabel('Error')
zlim([zMinError zMaxError])
set(gca,'fontsize',16,'Position',[0.55 0.03 0.4 0.30],'TitleFontWeight','normal','TitleFontSizeMultiplier',1.0)
box on

%%
subplot(3,2,3)
plot4=surf(zSave1)
zlabel('State')
zlim([zMin zMax])
set(gca,'fontsize',16,'Position',[0.08 0.35 0.4 0.30],'TitleFontWeight','normal','TitleFontSizeMultiplier',1.0)
box on

subplot(3,2,4)
plot1=surf(zSave1-zFull)
zlabel('Error')
zlim([zMinError zMaxError])
set(gca,'fontsize',16,'Position',[0.55 0.35 0.4 0.30],'TitleFontWeight','normal','TitleFontSizeMultiplier',1.0)
box on

%%
subplot(3,2,1)
surf(zCoarse)
zlabel('State')
zlim([zMin zMax])
set(gca,'fontsize',16,'Position',[0.08 0.67 0.4 0.30],'TitleFontWeight','normal','TitleFontSizeMultiplier',1.0)
box on

subplot(3,2,2)
surf(zCoarse-zFull)
zlabel('Error')
zlim([zMinError zMaxError])
set(gca,'fontsize',16,'Position',[0.55 0.67 0.4 0.30],'TitleFontWeight','normal','TitleFontSizeMultiplier',1.0)
box on

fig=gcf;
fig.PaperPosition(3:4)=figureSize;
saveas(gcf,'figure/figure.eps','epsc')
%%

figure(2)
clf
hold on

subplot('Position',[0.62 0.54 0.33 0.42])

surf(dWaveFull)
zlabel('$d_1 (x,y)$','Interpreter','latex')
zlim([dMin dMax])
set(gca,'fontsize',16,'TitleFontWeight','normal','TitleFontSizeMultiplier',1.0)
grid on
box on
subplot('Position',[0.62 0.08 0.33 0.42])

surf(dPeakFull)
zlabel('$d_2 (x,y)$','Interpreter','latex')
zlim([dMin dMax])
set(gca,'fontsize',16,'TitleFontWeight','normal','TitleFontSizeMultiplier',1.0)
grid on
box on

subplot('Position',[0.11 0.08 0.42 0.88])
surf(dFull)
zlabel('$d(x,y)=d_1 (x,y)+ d_2 (x,y)$','Interpreter','latex')
zlim([dMin dMax])
set(gca,'fontsize',16,'TitleFontWeight','normal','TitleFontSizeMultiplier',1.0)
grid on
box on

fig=gcf;
fig.PaperPosition(3:4)=[14 7];
saveas(gcf,'figure/spatial_figure2.eps','epsc')