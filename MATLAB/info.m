%% load data
clear all
clc
addpath(pwd);

addpath([pwd,'\matlab\mc']);
addpath([pwd,'\matlab\lattice2d']);

data_path = ['E:/visualize/inputInfo'];
% pathRegular = ['E:/visualize/inputHistogram/regular/'];
date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
save_address = ['E:/visualize/output/', date_time_string , '_Info'];
[status, msg, msgID] = mkdir(save_address)


%%---------------------------------------------------------------------------------
%% read data
Temperature = '0.849000';
DistancesSorted = h5read([data_path, '/NeighborsDistancesSorted.h5'],['/R0j']);
AcceptedProb = h5read([data_path, '/SelectedBondsProbability.h5'],['/', Temperature, '/AcceptedBondsProbability']);
SelectedProb = h5read([data_path, '/SelectedBondsProbability.h5'],['/', Temperature, '/SelectedBondsProbability']);

%%---------------------------------------------------------------------------------
%% find unique
[UniqueDistances, ~, idx] = uniquetol(DistancesSorted);
DistancesAccept = zeros(length(UniqueDistances), 1);
DistancesSelect = zeros(length(UniqueDistances), 1);

for i = 1:length(DistancesSorted)
    DistancesAccept(idx(i)) = DistancesAccept(idx(i)) + AcceptedProb(i);
    DistancesSelect(idx(i)) = DistancesSelect(idx(i)) + SelectedProb(i);
end

% xrange = 1:length(AcceptedProb);
cdfAccept = cumsum(AcceptedProb);
cdfSelect = cumsum(SelectedProb);
cdfAcceptDistance = cumsum(DistancesAccept);
cdfSelectDistance = cumsum(DistancesSelect);
%%---------------------------------------------------------------------------------
%% plot 
plot(UniqueDistances, cdfSelectDistance, 'Marker', '.', 'MarkerSize', 15);
% plot(UniqueDistances, cdfAcceptDistance, 'Marker', '.', 'MarkerSize', 15);
% xrange = 1:length(AcceptedProb);
xrange = 1:400;
% plot(xrange, AcceptedProb(xrange), 'Marker', '.', 'MarkerSize', 15);
% plot(xrange, SelectedProb, 'Marker', '.', 'MarkerSize', 15);
xlabel '$r$'
ylabel 'Number of resampled within distance'
% legendStrings = "$T$ = " + string(Temperature);
% legend(legendStrings, 'Location','northeast', ...
%     'Interpreter','latex')
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
% cmap = linspecer(length(Temperature));
% ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'resampled.png']);

%%---------------------------------------------------------------------------------
%% plot cdf
cdfAccept
% plot(UniqueDistances, DistancesSelect, 'Marker', '.', 'MarkerSize', 15);
% plot(UniqueDistances, DistancesAccept, 'Marker', '.', 'MarkerSize', 15);
% xrange = 1:length(AcceptedProb);
xrange = 1:400;
plot(xrange, AcceptedProb(xrange), 'Marker', '.', 'MarkerSize', 15);
% plot(xrange, SelectedProb, 'Marker', '.', 'MarkerSize', 15);
xlabel '$r$'
ylabel 'Number of Accepted'
% legendStrings = "$T$ = " + string(Temperature);
% legend(legendStrings, 'Location','northeast', ...
%     'Interpreter','latex')
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
% cmap = linspecer(length(Temperature));
% ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'm.png']);

