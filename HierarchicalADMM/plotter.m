set(0,'defaultAxesFontSize',11,...
      'defaultTextInterpreter','latex', ...
      'defaultLegendInterpreter', 'latex',...
      'defaultLineLineWidth', 1.0,...
      'defaultFigureVisible','off')
yl={'Primal Residual','Dual Residual',...
    'Objective','Aug. Lagrangian',...
    'Lyapunov Function'};

tbl = [];
tbl2= [];

fg=figure();
fg.Position(3:4)=[1200 550];
for i=2:6
    tbl_dec = csvread(sprintf('output/dec-%i.csv',i));
    tbl_mul = csvread(sprintf('output/mul-%i.csv',i));
    stat = csvread(sprintf('output/stat-%i.csv',i));

    % extra = [stat(8) stat(9) stat(12) stat(12)];
    extra = [max(stat(8),stat(10)) max(stat(9),stat(11)) stat(12) stat(12)];
    
    l_dec=size(tbl_dec,1)-1;
    iter_dec = 0:l_dec;
    l_mul=size(tbl_mul,1)-1;
    iter_mul = 0:l_mul;

    % wrt iteration steps
    for j=1:4
        sp=subplot(4,5,i-1+5*(j-1));
        hold on; box on; grid on; set(gca,'YScale','log');
        plot(0:max(l_dec,l_mul),ones(max(l_dec,l_mul)+1,1)*extra(j),'k--');
        plot(iter_mul,tbl_mul(:,j),'r-','LineWidth',2);
        plot(iter_dec,tbl_dec(:,j),'b-');
        xlabel('Iteartion Step');
        ylabel(yl{j});
        xlim([0 max(l_dec,l_mul)]);
        if j >= 3
            if i==2
                dd=.01;
                uu=.02;
            else
                dd=.2;
                uu=.4;
            end
            sp.YLim = [10^(log10(extra(j))-dd) 10^(log10(extra(j))+uu)];
        end
    end
    % saveas(fg,sprintf('figure/figure%i.eps',i),'epsc');

    % wrt compuqtion time
    % fg=figure();
    % fg.Position(3:4)=[1200 250];
    % for j=1:4
    %     subplot(1,4,j);
    %     hold on; box on; grid on; set(gca,'YScale','log');
    %     plot([0 tbl_dec(end,6)],ones(2,1)*extra(j),'k--');
    %     plot(tbl_mul(:,6)+stat(2),tbl_mul(:,j),'r-','LineWidth',2);
    %     plot(tbl_dec(:,6),tbl_dec(:,j),'b-');
    %     xlabel('Time (sec)');
    %     ylabel(yl{j});
    %     xlim([0 max(tbl_dec(end,6),tbl_mul(end,6))]);
    % end
    % saveas(fg,sprintf('figure/figure%i-t.eps',i),'epsc');

    tbl = [tbl ; stat(4:7)'];
    tbl2= [tbl2; [stat(1) stat(12)/1e6 ...
                  tbl_dec(end,end) l_dec tbl_dec(end,3)/1e6 ...
                  stat(2) tbl_mul(end,end)+stat(2) l_mul tbl_mul(end,3)/1e6 ...
                 stat(3)]];
end