% Sungho Shin (sungho.shin@wisc.edu) and Victor Zavala (zavalatejeda@wisc.edu)
% Makes posterior plots
set(0,'defaultFigureVisible','off')

set(0, 'DefaultFigureRenderer', 'painters');

close all
clf
clear
clc

distmax=2.0;
nSpecies=12;
nMesh=40;

numOfParams=nSpecies*(nSpecies+1);

outputPath=['../output/UQ/posterior/'];
figurePath='../figure/';

mkdir(figurePath);
fg=figure;
fg.Position(3:4)=[1920 984];
clf

map=hsv;
for i=1:nSpecies
    clrs{i}=map(floor(i*(64/nSpecies)),:)*0.6;
end
greyclr=[.5 .5 .5];
params=[];
while i<5000
    i=i+1;
    filePath=[outputPath 'param' num2str(i) '.csv'];
    if exist(filePath)~=0
        params=[params csvread(filePath)];
    end
end

covar=cov(params');
meanPar=mean(params,2);
varPar=var(params')';
std=sqrt(varPar);
rsd=abs(sqrt(varPar)./meanPar);


for i=1:nSpecies
    for j=1:nSpecies+1
        k=(nSpecies+1)*(i-1)+j;
        distmax = 3*std(k);
        subplot(nSpecies,nSpecies+1,13*(i-1)+j);
        h=histogram(params(k,:),linspace(meanPar(k)-distmax,meanPar(k)+distmax,50));
        h.EdgeColor='none';
        h.FaceAlpha=1;
        h.FaceColor=clrs{i};
        h.Parent.YAxis.TickLabels={};
        xlim([meanPar(k)-distmax meanPar(k)+distmax]);
        ylim([0 50])
        h.Parent.XAxis.TickLabels={};
    end
end

speciesOrder={'BH','CA','BU','PC','BO','BV','BT','EL','FP','CH','DP','ER'};
for i=1:nSpecies
    vtext = uicontrol('style','text');
    set(vtext,'String',speciesOrder{i});
    vtext.Units='normalized';
    vtext.Position(1:2)=[0.105,0.058+(13-i)*(0.835/nSpecies)];
    vtext.Position(3)=vtext.Position(3)*0.8;
    vtext.Position(4)=vtext.Position(4)*1.3;
    vtext.FontSize=24;
    vtext.ForegroundColor=clrs{i};
    vtext.BackgroundColor=[1 1 1];
    
    htext = uicontrol('style','text');
    set(htext,'String',speciesOrder{i});
    htext.Units='normalized';
    htext.Position(1:2)=[0.081+(i+1)*(0.788/(nSpecies+1)),0.927];
    htext.Position(3)=htext.Position(3)*0.8;
    htext.Position(4)=vtext.Position(4)*1.15;
    htext.FontSize=24;
    htext.ForegroundColor=[0 0 0];
    htext.BackgroundColor=[1 1 1];
end

saveas(fg,[figurePath 'UQ_conf.eps'],'epsc')

fg=figure;
fg.Position(3:4)=[1920 984];
clf


for ii=1:nSpecies
    for jj=1:nSpecies
        i=(nSpecies+1)*(ii-1)+1;
        j=(nSpecies+1)*(ii-1)+1+jj;


        sp=subplot(nSpecies,nSpecies,nSpecies*(ii-1)+jj);
        A=inv([covar(i,i) covar(i,j); covar(j,i) covar(j,j)]);
        m=[meanPar(i);meanPar(j)];
        b=5.991;
        invA=inv(A);
        xlim1=-sqrt(b*invA(1,1))+m(1);
        xlim2=+sqrt(b*invA(1,1))+m(1);
        ylim1=-sqrt(b*invA(2,2))+m(2);
        ylim2=+sqrt(b*invA(2,2))+m(2);
        
        x1=linspace(xlim1,xlim2,nMesh);
        x2=linspace(ylim1,ylim2,nMesh);
        [X1,X2]=meshgrid(x1,x2);
        Z=zeros(nMesh);
        for k=1:nMesh
            for l=1:nMesh
                X=[X1(k,l);X2(k,l)];
                Z(k,l)=(X-m)'*A*(X-m);
            end
        end
        hold on
        plot(params(i,:),params(j,:),'.','Color',clrs{ii},'MarkerSize',5) ;
        ctr=contour(X1,X2,Z,[b,b],'Color','k','LineWidth',1.5);
        sp.XTick=[];
        sp.YTick=[];
        
        box on
    end
end

labelV={'BH','CA','BU','PC','BO','BV','BT','EL','FP','CH','DP','ER'};
labelH={'BH','CA','BU','PC','BO','BV','BT','EL','FP','CH','DP','ER'};
for i=1:nSpecies
    vtext = uicontrol('style','text');
    set(vtext,'String',labelV{i});
    vtext.Units='normalized';
    vtext.Position(1:2)=[0.105,0.055+(13-i)*(0.835/nSpecies)];
    vtext.Position(3)=vtext.Position(3)*0.65;
    vtext.Position(4)=vtext.Position(4)*1.3;
    vtext.FontSize=24;
    vtext.ForegroundColor=clrs{i};
    vtext.BackgroundColor=[1 1 1];
    
    htext = uicontrol('style','text');
    set(htext,'String',labelH{i});
    htext.Units='normalized';
    htext.Position(1:2)=[0.075+i*(0.792/nSpecies),0.93,];
    htext.Position(3)=htext.Position(3)*0.8;
    htext.Position(4)=htext.Position(4)*1.3;
    htext.FontSize=24;
    htext.ForegroundColor=[0 0 0];
    htext.BackgroundColor=[1 1 1];
end

saveas(fg,[figurePath 'UQ_ellipse.eps'],'epsc')

set(0, 'defaultTextInterpreter','latex','defaultLegendInterpreter', ...
       'latex','defaultAxesTickLabelInterpreter','latex')

fg= figure;
fg.Position(3:4)=[600 400];
std_rsh=reshape(std,13,12)';
h=heatmap(std_rsh');
h.YData = {'-',speciesOrder{:}};
h.XData = speciesOrder;
h.ColorLimits = [0 max(abs(std_rsh(:)))];
h.Colormap=[ ones(64,1) linspace(1,0,64)' linspace(1,0,64)'];
h.FontSize=20;
h.CellLabelColor='none';
saveas(fg,[figurePath 'UQ_hmap_std.eps'],'epsc')

for i=1:156
    ec{i} = ' ';
end

fg= figure;
cmap =[ ones(32,1) linspace(0,1,32)' linspace(0,1,32)';
        linspace(1,0,32)' linspace(1,0,32)' ones(32,1)];
fg.Position(3:4)=[600 400];
pcor = covar./(std*std');
h=heatmap(pcor);
h.ColorLimits = [-max(abs(pcor(:))) max(abs(pcor(:)))];
h.Colormap=cmap;
h.FontSize=20;
h.CellLabelColor='none';
h.GridVisible='off';
h.YDisplayLabels=ec;
h.XDisplayLabels=ec;
saveas(fg,[figurePath 'UQ_hmap_pcor.eps'],'epsc')

% Third momdent
mom3 = moment(params',3)'./(std.^3);
mom3_rsh=reshape(mom3,13,12)';
cmap =[ ones(32,1) linspace(0,1,32)' linspace(0,1,32)';
        linspace(1,0,32)' linspace(1,0,32)' ones(32,1)];

fg= figure;
fg.Position(3:4)=[600 400];
h=heatmap(mom3_rsh');
h.YData = {'-',speciesOrder{:}};
h.XData = speciesOrder;
h.ColorLimits = [-max(abs(mom3_rsh(:))) max(abs(mom3_rsh(:)))];
h.Colormap=cmap;
h.FontSize=20;
h.CellLabelColor='none';
saveas(fg,[figurePath 'UQ_hmap_mom3.eps'],'epsc')

% Forth momdent
mom4 = moment(params',4)'./(std.^4);
mom4_rsh=reshape(mom4,13,12)';
cmap =[ ones(32,1) linspace(0,1,32)' linspace(0,1,32)';
        linspace(1,0,32)' linspace(1,0,32)' ones(32,1)];

fg= figure;
fg.Position(3:4)=[600 400];
h=heatmap(mom4_rsh');
h.YData = {'-',speciesOrder{:}};
h.XData = speciesOrder;
h.ColorLimits = [3-max(abs(3-mom4_rsh(:))) 3+max(abs(3-mom4_rsh(:)))];
h.Colormap=cmap;
% h.Colormap=[ ones(64,1) linspace(1,0,64)' linspace(1,0,64)'];
h.FontSize=20;
h.CellLabelColor='none';
saveas(fg,[figurePath 'UQ_hmap_mom4.eps'],'epsc')