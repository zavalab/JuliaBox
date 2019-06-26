data = load('daily_prediction.csv');
day = data(:,1);
X   = data(:,2);
P   = data(:,3);
TSI_X = 10*(6-(2.04-0.68*log(X))/log(2));
TSI_P = 10*(6-log(48./P)/log(2));

figure()
plot([0;day],50*ones(1+length(day),1),'-.k','LineWidth',1);
hold on

plot([0;day],10*ones(1+length(day),1),'-.k','LineWidth',1);
plot([0;day],0*ones(1+length(day),1),'-.k','LineWidth',1);
plot([0;day],[X(1);X],'k','LineWidth',2);
axis([0 210 0 60]);
%ylabel('Concentration of Chl-a (mg/m^3)');
x = 3;
y = 56;
% txt = '(Scenario I)';
% text(x,y,txt,'HorizontalAlignment','left')
XAxisLine = 'off';
x0=10;
y0=10;

width=600;
height=337.5;
set(gcf,'units','points','position',[x0,y0,width,height]);
set(gca,'xtick',[]);
xt = get(gca, 'YTick');
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