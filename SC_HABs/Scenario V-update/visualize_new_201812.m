data = load('daily_chla_distribution.csv');
%data2 = load('daily_prediction.csv');
day = data(2:end,1);
chla = data(2:end,2:end);
%X   = data2(:,2);
%P   = data2(:,3);
%TSI_X = 10*(6-(2.04-0.68*log(X))/log(2));
%TSI_P = 10*(6-log(48./P)/log(2));

chla_surface = 1/5*(chla(:,1) + chla(:,2) + chla(:,3) + chla(:,4) + chla(:,5));

figure()
plot([0;day],50*ones(1+length(day),1),'r--','LineWidth',2);
hold on

plot([0;day],10*ones(1+length(day),1),'--','Color',[0.95 0.7 0.05],'LineWidth',2);
plot([0;day],0*ones(1+length(day),1),'g--','LineWidth',2);
plot([0;day],[chla_surface(1);chla_surface],'b','LineWidth',2);
axis([0 210 0 150]);
%ylabel('Concentration of Chl-a (mg/m^3)');
x = 3;
y = 56;
% txt = '(Scenario I)';
% text(x,y,txt,'HorizontalAlignment','left')
XAxisLine = 'off';
x0=10;
y0=10;

width=600;
height=200;
set(gcf,'units','points','position',[x0,y0,width,height]);
set(gca,'xtick',[]);
xt = get(gca, 'YTick');
set(gca, 'FontSize', 14);
hold off

figure()
colormap('jet');
imagesc(1:210,0:0.1:2.5,(chla(:,1:25))');
hold on
%colorbar('Ticks',[0,50,100,150,200,250,300]);
caxis([0 300])
width=600;
height=200;
set(gcf,'units','points','position',[x0,y0,width,height]);
set(gca,'xtick',[]);
%set(gca,'ytick',[]);
%xt = get(gca, 'YTick');
%yt = get(gca, 'XTick');
set(gca, 'FontSize', 14);

hold off

% figure()
% plot(92:183,TSI_X(92:183),'g','LineWidth',2);
% hold on
% plot(92:183,TSI_P(92:183),'r','LineWidth',2);
% scatter([92+7;92+9;92+18;92+23;123+5;123+6;123+19;123+22;123+30;154+5;154+24],[53;58;54;58;57;56;57;55;55;54;54],'b','LineWidth',1.5);
% axis([92 183 30 80]);
% x = 93.5;
% y = 77;
% txt = '(Scenario I)';
% text(x,y,txt,'HorizontalAlignment','left')
% XAxisLine = 'off';
% x0=10;
% y0=10;
% 
% width=600;
% height=200;
% set(gcf,'units','points','position',[x0,y0,width,height]);
% set(gca,'xtick',[]);
% xt = get(gca, 'YTick');
% set(gca, 'FontSize', 14);
% hold off