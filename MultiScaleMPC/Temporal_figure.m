close all
mkdir figure


clm=colormap();
clm2=colormap('hsv');

timeGrid=csvread('output/Temporal/timeGrid.csv');
dFull=csvread('output/Temporal/dFull.csv');
d1=csvread('output/Temporal/d1.csv');
d2=csvread('output/Temporal/d2.csv');
zFull=csvread('output/Temporal/zFull.csv');
zGS=csvread('output/Temporal/zGS.csv');
zCoarse=csvread('output/Temporal/zCoarse.csv');
zSave1=csvread('output/Temporal/zSave1.csv');
zSave2=csvread('output/Temporal/zSave10.csv');
zSaveCon=csvread('output/Temporal/zSave30.csv');
errorSave=csvread('output/Temporal/errorSave.csv');
errorSaveWC=csvread('output/Temporal/errorSaveWC.csv');
zMax=max([max(zFull(:)),max(zSave1(:)),max(zSave2(:)),max(zCoarse(:))]);
zMin=min([min(zFull(:)),min(zSave1(:)),min(zSave2(:)),min(zCoarse(:))]);
dMax=max(dFull(:))+2;
dMin=min(dFull(:))-2;

figure(1)
clf
hold on
sp1=subplot(1,2,1)
plot1=plot(timeGrid,dFull,'Color',clm(12,:),'LineWidth',1.5);
xlabel('Time')
lab=ylabel('$d(t)=d_1(t)+d_2(t)$','Interpreter','Latex');
ylim([dMin dMax])
set(gca,'fontsize',14)
grid on
box on

sp2=subplot(2,2,2)
plot1=plot(timeGrid,d1,'Color',clm(12,:),'LineWidth',1.5);
xlabel('Time')
ylim([dMin dMax])
ylabel('$d_1(t)$','Interpreter','Latex')
set(gca,'fontsize',14)
grid on
box on

sp4=subplot(2,2,4)
plot1=plot(timeGrid,d2,'Color',clm(12,:),'LineWidth',1.5);
xlabel('Time')
ylim([dMin dMax])
ylabel('$d_2(t)$','Interpreter','Latex')
set(gca,'fontsize',14)
grid on
box on

sp1.Position=[.09 .09 .42 .86];
sp2.Position=[.57 .57 .42 .38];
sp4.Position=[.57 .09 .42 .38];

fig=gcf;
fig.PaperPosition(3:4)=[14 6];
saveas(gcf,'figure/temporal_figure1.eps','epsc')

figure(2)
clf
sp1=subplot(2,2,1);

hold on
grid on
plot(timeGrid,zCoarse,'d-','Color',clm(12,:),'LineWidth',1.5)
plot(timeGrid,zFull,'Color','k','LineWidth',1.5)

xlabel('Time')
ylabel('State')
ylim([zMin zMax])
leg1=legend('Coarse','Optimal');
leg1.Position(1:2)=[.08 .58];
set(sp1,'fontsize',14)
box on

sp2=subplot(2,2,2);
hold on
grid on
plot(timeGrid,zSave1,'d-','Color',clm(12,:),'LineWidth',1.5)
plot(timeGrid,zFull,'k','LineWidth',1.5)

xlabel('Time')
ylabel('State')
ylim([zMin zMax])
leg2=legend('GS #1','Optimal');
leg2.Position(1:2)=[.56 .58];
set(sp2,'fontsize',14)
box on

sp3=subplot(2,2,3);
hold on
grid on
plot(timeGrid,zSaveCon,'d-','Color',clm(12,:),'LineWidth',1.5)
plot(timeGrid,zFull,'k','LineWidth',1.5)

xlabel('Time')
ylabel('State')
ylim([zMin zMax])
leg3=legend('GS #30','Optimal');
leg3.Position(1:2)=[.08 .1];
set(sp3,'fontsize',14)
box on

% sp4=subplot(2,2,4);
% hold on
% grid on
% plot(timeGrid,zSave1noCoarse,'d-','Color',clm2(2,:),'LineWidth',1.5)
% plot(timeGrid,zFull,'k','LineWidth',1.5)

% xlabel('Time')
% ylabel('State')
% ylim([zMin zMax])
% leg3=legend('GS #1(No Coarse)','Optimal');
% leg3.Position(1:2)=[.549 .1];
% set(sp3,'fontsize',14)
% box on
% set(sp4,'fontsize',14)

fig=gcf;
fig.PaperPosition(3:4)=[14 12];

sp1.Position=[.09 .57 .4 .4];
sp2.Position=[.57 .57 .4 .4];
sp3.Position=[.09 .09 .4 .4];
sp4.Position=[.57 .09 .4 .4];

saveas(gcf,'figure/temporal_figure2.eps','epsc')

