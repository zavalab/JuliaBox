figurePath='../figure/';
nSpecies=12;

mkdir(figurePath);
fg=figure(1);
fg.Position(3:4)=[1200 200];
clf;hold on;

map=hsv;
for i=1:nSpecies
    clrs{i}=map(floor(i*(64/nSpecies)),:)*0.6;
end
greyclr=[.5 .5 .5];
markers={'-',':'};

% Saturable
folderNames={'saturable','standard'};
exps={'M','M','P','P'};
spind={1,1,1,2};
inds={1,9,4,60};
hlind={[],[],[],[]};

% CVaR
% folderNames={'cvar','standard'};
% exps={'P','P','P','P'};
% spind={1,2,2,2};
% inds={4,60,14,16};
% hlind={5,2,6,2};

for k=1:4
    subplot(1,4,k);set(gca, 'FontSize', 14)
    hold on; box on; grid on;

    exp=exps{k};
    outputNumber=inds{k};
    j=spind{k};

    outputPath=['../output/' folderNames{1} '/yout/' exp '/y' num2str(outputNumber) '/'];

    species=csvread([outputPath 'species.csv']);
    timePoints=csvread([outputPath 'timePoints.csv']);
    timePointsExtended=csvread([outputPath 'timePointsExtended.csv']);
    numSpecies=length(species);

    yexp=csvread([outputPath 'yexp' num2str(j) '.csv']);
    plot(timePoints,yexp,'*','Color',greyclr)
    plot(timePoints(hlind{k}),yexp(hlind{k}),'ro','MarkerSize',12)
    for kk=1:2
        folderName=folderNames{kk};
        outputPath=['../output/' folderName '/yout/' exp '/y' num2str(outputNumber) '/'];
        yout=csvread([outputPath 'yout' num2str(j) '.csv']);
        plot(timePointsExtended,yout,'Color', ...
             clrs{species(j)},'LineWidth',2,'LineStyle',markers{kk})
    end
    xlim([0 timePoints(end)])
end
saveas(figure(1),[figurePath 'fit_select.eps'],'epsc')