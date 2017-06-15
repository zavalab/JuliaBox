clc
clear all;
close all;

% Reading values for tradeoff
fidi = fopen('tradeoff.csv');
d = textscan(fidi, '%f %f %f %f %f %f %f %f', 'Delimiter',',', 'HeaderLines',3);
fclose(fidi);

epsilon = cell2mat(d(1));;
invcost = cell2mat(d(2));;
transcost = cell2mat(d(3));
opcost = cell2mat(d(4));
daily_cost = cell2mat(d(5));
umanure = cell2mat(d(6));
umanure_percent = cell2mat(d(7));

%scatter(transcost(8:end), umanure_percent(8:end))
%hold on;
%xlabel('Transportation Cost ($/day)');
%ylabel('Percentage of Manure Unprocessed (kg/day)');
%legend('Budget: $ 5 million')


f = fit(daily_cost, umanure_percent, 'linearinterp');
plot(f);
hold on;
%scatter(daily_cost, umanure_percent);
xlabel('d_{cost} ($/day)');
ylabel('\phi (%)');
%grid on;
legend('Pareto Front ');