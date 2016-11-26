% plot compromise solution of Formulation D (nomalized) + Utopian points
% Yankai Cao, Siyu Chen, Luis Fuentes
% UW-Madison, 2016
Costref = 1;
GHGEref = 1;
SWref = 1;

%read data, the original data is normalized
Cost_UP = dlmread('Cost_UP.txt')*Costref;
GHGE_UP = dlmread('GHGE_UP.txt')*GHGEref;
SW_UP = dlmread('SW_UP.txt')*SWref;
Cost = dlmread('Cost_ao0.999_as0.txt')*Costref;
GHGE = dlmread('GHGE_ao0.999_as0.txt')*GHGEref;
SW = dlmread('SW_ao0.999_as0.txt')*SWref;

h = figure;
set(gca,'fontsize',18)
scatter3( Cost_UP, GHGE_UP, SW_UP, 'blue');
hold on
set(gca,'fontsize',18)
scatter3( Cost, GHGE, SW, 'red');
grid on
axis([1 2.6 1 1.9  1 2.1])

%add legend and labels
l = legend('UP', 'Minimize $\langle\!\langle \mathbf{F}(x) \rangle\!\rangle_{\alpha_{\mathcal{S}}=0,\alpha_\mathcal{O}=1}$','Location','northeast');
set(l,'Interpreter','latex') 
xlabel('Normalized Cost')
ylabel('Normalized GHGE Emission')
zlabel('Normalized Water demand')
