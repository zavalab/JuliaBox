% Sungho Shin (sungho.shin@wisc.edu) and Victor Zavala (zavalatejeda@wisc.edu)
% Makes scalability plots (need to enter the data manually)

clear
clc
sol_time_PIPS=[6982.059789177 3868.111628927 2081.580929735 1125.269504376 596.146505297]/60;
sol_time_IPOPT=1769.632958108/60;
cores=[1 2 4 8 16];

fg=figure(1)
clf
subplot(1,2,1)
set(gca, 'FontSize', 14)
hold on
plot(cores,sol_time_IPOPT*ones(length(sol_time_PIPS),1),'b-')
plot(cores,sol_time_PIPS,'dr-')
lg=legend('IPOPT','PIPS-NLP');
lg.FontSize=14;
grid on
box on
xlabel('Number of Processors','FontSize',14)
ylabel('Solution Time (min)','FontSize',14)
xlim([1 16])

subplot(1,2,2)
hold on
plot(cores,sol_time_PIPS(1)./sol_time_IPOPT*ones(length(sol_time_PIPS),1),'b-')
plot(cores,sol_time_PIPS(1)./sol_time_PIPS,'dr-')
pl=plot(cores,cores,'Color',[.5 .5 .5],'LineStyle','--')
lg=legend('','');
lg.Location='northwest';
lg.FontSize=14;
lg.String={'IPOPT','PIPS-NLP'}
ylim([1 16])
xlim([1 16])

grid on
box on
xlabel('Number of Processors','FontSize',14)
ylabel('Speed Up','FontSize',14)

fg.Position(3:4)=[600 250];

set(gca, 'FontSize', 14)
saveas(fg,'../figure/scl_P2.eps','epsc')

fg=figure(2)
clf
comm_size=[12,24,36,48];
n_var=[ 83004,340824,773460,1380912];
sol_time=[31.794692168,128.802693354,447.63088772,868.169078085]/60;

clf
subplot(1,2,1)
set(gca, 'FontSize', 14)
hold on
pl=plot(comm_size,n_var,'kd-')
lg.FontSize=14;
grid on
box on
xlabel('Community Size','FontSize',14)
ylabel('Number of Variables','FontSize',14)

subplot(1,2,2)
hold on
plot(comm_size,sol_time,'kd-')
lg.Location='northwest';
lg.FontSize=14;

grid on
box on
xlabel('Community Size','FontSize',14)
ylabel('Solution Time (min)','FontSize',14)

fg.Position(3:4)=[600 250];

set(gca, 'FontSize', 14)
saveas(fg,'../figure/scl_S.eps','epsc')
