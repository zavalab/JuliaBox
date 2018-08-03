% Sungho Shin (sungho.shin@wisc.edu) and Victor Zavala (zavalatejeda@wisc.edu)
% Makes error histogram

set(0, 'DefaultFigureRenderer', 'painters');
set(0, 'defaultTextInterpreter','latex','defaultLegendInterpreter', ...
       'latex','defaultAxesTickLabelInterpreter','latex')


close all
clear
clc
clf

quantile{1}=0;
quantile{2}=0.9;
quantile{3}=0.9;
quantile{4}=0.9;

mses_{1}=sort(csvread('../output/standard/mses.csv'))/2;
mses_{2}=sort(csvread('../output/cvar/mses.csv'))/2;

k=floor(length(mses_{1})/2);
xmin=1e-4/2;
xmax=max([mses_{1};mses_{2}]);
xmin2=1e2/2;
ymax=200;
ymax2=20;

edges=logspace(log10(xmin),log10(xmax),100);
edges2=logspace(log10(xmin2),log10(xmax),40);


for i=1:2
    fg=figure;
    histogram(mses_{i},edges,'FaceColor','black');
    hold on
    
    m1{i}=mean(mses_{i});
    m2{i}=mean(mses_{i}(ceil(quantile{2}*length(mses_{i})):end))
    m3{i}=mean(mses_{i}(ceil(quantile{3}*length(mses_{i})):end));
    m4{i}=mean(mses_{i}(ceil(quantile{4}*length(mses_{i})):end));
    
    plot(ones(1,20)*m1{i},linspace(0,ymax,20),'r','LineWidth',4);
    plot(ones(1,20)*m4{i},linspace(0,ymax,20),'b','LineWidth',4);
    
    % xlabel('Error','FontSize',13);
    lg=legend('Errors','Mean','CVaR$_{0.9}$');
    lg.FontSize=11;
    lg.Location='northwest';
    xlim([xmin xmax]);
    ylim([0 ymax]);
    grid on
    box on
    set(gca,'xscale','log')
    fg.Position(3:4)=[400 180];
    saveas(fg,['../figure/mse_histogram' num2str(i) '.eps'],'epsc')

    fg=figure;
    mse_zoom=mses_{i}(mses_{i}>=xmin*(xmax/xmin)^.8);
    histogram(mse_zoom,edges2,'FaceColor','black');
    hold on
    
    plot(ones(1,20)*m4{i},linspace(0,ymax2,20),'b','LineWidth',4);
    
    % xlabel('Error','FontSize',13);
    lg=legend('Errors','CVaR$_{0.9}$');
    lg.FontSize=11;
    lg.Location='northeast';
    xlim([xmin2 xmax]);
    ylim([0 ymax2]);
    grid on
    box on
    set(gca,'xscale','log')
    fg.Position(3:4)=[400 180];
    saveas(fg,['../figure/mse_histogram_zoom' num2str(i) '.eps'],'epsc')
end

