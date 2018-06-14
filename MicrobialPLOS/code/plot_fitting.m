% Sungho Shin (sungho.shin@wisc.edu) and Victor Zavala (zavalatejeda@wisc.edu)
% Makes model fitting plots

ind=1;
nSpecies=12;
folderName='standard';
figurePath='../figure/';
mkdir(figurePath);
fg=figure(1);
ymax=1.5;
clf
hold on

map=hsv;
for i=1:nSpecies
    clrs{i}=map(floor(i*(64/nSpecies)),:)*0.6;
end
greyclr=[.5 .5 .5];

for monoind=1:nSpecies
    exp='M';
    outputNumber=monoind;

    outputPath=['../output/' folderName '/yout/' exp '/y' num2str(outputNumber) '/'];

    species=csvread([outputPath 'species.csv']);
    timePoints=csvread([outputPath 'timePoints.csv']);
    timePointsExtended=csvread([outputPath 'timePointsExtended.csv']);
    numSpecies=length(species);

    yexp=zeros(numSpecies,length(timePoints));
    yout=zeros(numSpecies,length(timePointsExtended));
    for j=1:numSpecies
        yexp(j,:)=csvread([outputPath 'yexp' num2str(j) '.csv']);
        yout(j,:)=csvread([outputPath 'yout' num2str(j) '.csv']);
    end
    yexp;

    for j=1:numSpecies
        subplot(nSpecies,nSpecies,nSpecies*(species-1)+species);
        plot(timePoints,yexp(j,:),'*','Color',greyclr);
        hold on
        plot(timePointsExtended,yout(j,:),'Color',clrs{species},'LineWidth',2);
        box on
        grid on
        xlim([0 timePoints(end)]);
    end
end

for pairid=1:66
    exp='P';
    outputNumber=pairid;

    outputPath=['../output/' folderName '/yout/' exp '/y' num2str(outputNumber) '/'];
    figurePath='../figure/';

    species=csvread([outputPath 'species.csv']);
    timePoints=csvread([outputPath 'timePoints.csv']);
    timePointsExtended=csvread([outputPath 'timePointsExtended.csv']);
    numSpecies=length(species);

    yexp=zeros(numSpecies,length(timePoints));
    yout=zeros(numSpecies,length(timePointsExtended));
    for j=1:numSpecies
        yexp(j,:)=csvread([outputPath 'yexp' num2str(j) '.csv']);
        yout(j,:)=csvread([outputPath 'yout' num2str(j) '.csv']);
    end
    yexp;

    for j=1:numSpecies
        if j==1
            subplot(nSpecies,nSpecies,nSpecies*(species(1)-1)+species(2))
        else
            subplot(nSpecies,nSpecies,nSpecies*(species(2)-1)+species(1))
        end
        plot(timePoints,yexp(j,:),'*','Color',greyclr)
        hold on
        plot(timePointsExtended,yout(j,:),'Color',clrs{species(j)},'LineWidth',2)
        grid on
        box on
        xlim([0 timePoints(end)])
    end
end

speciesOrder={'BH','CA','BU','PC','BO','BV','BT','EL','FP','CH','DP','ER'};
for i=1:nSpecies
    vtext = uicontrol('style','text');
    set(vtext,'String',speciesOrder{i});
    vtext.Units='normalized';
    vtext.Position(1:2)=[0.098,0.06+(13-i)*(0.835/nSpecies)];
    vtext.Position(3)=vtext.Position(3)*0.75;
    vtext.Position(4)=vtext.Position(4)*1.3;
    vtext.FontSize=24;
    vtext.ForegroundColor=clrs{i};
    vtext.BackgroundColor=[1 1 1];
    
    htext = uicontrol('style','text');
    set(htext,'String',speciesOrder{i});
    htext.Units='normalized';
    htext.Position(1:2)=[0.08+i*(0.792/nSpecies),0.93,];
    htext.Position(3)=htext.Position(3)*0.75;
    htext.Position(4)=htext.Position(4)*1.3;
    htext.FontSize=24;
    htext.ForegroundColor=[0 0 0];
    htext.BackgroundColor=[1 1 1];
end

fg.Position(3:4)=[1920 984];
saveas(fg,[figurePath 'fit2.eps'],'epsc')