%% load data
clear all
clc
addpath(pwd);

addpath([pwd,'\matlab\mc']);
addpath([pwd,'\matlab\lattice2d']);

path = ['E:/visualize/inputNonEquilibrium/'];
% pathRegular = ['E:/visualize/inputHistogram/regular/'];
date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
save_address = ['E:/visualize/output/', date_time_string , '_NonEquilibrium'];
[status, msg, msgID] = mkdir(save_address)

% Llist = readmatrix([path, 'list.txt'])
% Nlist = Llist.^2;
% N = length(Nlist);

% Templist = readmatrix([path, 'Templist.txt'])
% NTemp = length(Templist);
% data = readmatrix([path 'MeanX.csv']);
% error = readmatrix([path 'errorX.csv']);
myxLabel = '$T$';
% for i=1:N
%     EnsembleIndex = num2str(i);
%     data_path = [pathRegular, EnsembleIndex , '/'];
% 
%     % [status, msg, msgID] = mkdir(folder_address)
% 
% %     dataRegular(:,:,i) = dlmread([data_path 'MeanX.txt']);
% %     errorRegular(:,:,i) = readmatrix([data_path 'errorX.txt']);
% end

dpi = 600;

% data = dlmread([path 'MeanX.txt']);
% error = dlmread([path 'errorX.txt']);

% orderParam = data(:,1);
% myxLabel = '$B\mu /\Lambda $';


% error = error / sqrt(10); % for results before 2022.08.25


% locations = readmatrix([path, '0/locations.txt']);
% DistanceVec = readmatrix([path, 'DistanceVec.txt']);
% % Plabel = readmatrix([path, 'Plabel.txt']);
% Mean_GVec = readmatrix([path, 'Mean_GVec.txt']);
% Mean_GConVec = readmatrix([path, 'Mean_GConVec.txt']);
% error_GVec = readmatrix([path, 'error_GVec.txt']);
% error_GConVec = readmatrix([path, 'error_GConVec.txt']);


% for i=1:Nensemble
%     EnsembleIndex = num2str(i-1);
%     data_path = [path, EnsembleIndex , '/'];
%     
%     EtimeSeries(:,i) = dlmread([data_path 'EtimeSeries.txt']);
%     mtimeSeries(:,i) = dlmread([data_path 'mtimeSeries.txt']);
%     mpPlanetimeSeries(:,i) = dlmread([data_path 'mPlanetimeSeries.txt']);
%     
% %     EtimeSeries(:,i) = readmatrix([data_path 'EtimeSeries.csv']);
% %     mtimeSeries(:,i) = readmatrix([data_path 'mtimeSeries.csv']);
% %     mpPlanetimeSeries(:,i) = readmatrix([data_path 'mPlanetimeSeries.csv']);
% end
% %%--------------------------------------------------------------------------------------
% %% save address (for loaded data)
% date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
% save_address = ['d:/visualize/output/', date_time_string , '_Histogram'];
% [status, msg, msgID] = mkdir(save_address)
%%--------------------------------------------------------------------------------------
%% Read Data
filelist = dir(fullfile(path, '**\results.h5'));  %get list of files named "results.h5" in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list

start = 1;  % start index
max_count = realmax;
for i = 1:length(filelist)
    data_path = [filelist(i).folder, '\', filelist(i).name];
    Temporary = h5read(data_path, '/Magnetization/mean');
    max_count = min(max_count, length(Temporary));
end

% max_count = double(max_count);
count = max_count
t_max = max_count
% count = 100;  % number of elements to read
Nbin = 5;
for i = 1:length(filelist)
    data_path = [filelist(i).folder, '\', filelist(i).name];
    % Read data from the HDF5 file
    meanMagnetization(:, i) = h5read(data_path, '/Magnetization/mean', start, count);
    errorMagnetization(:, i) = h5read(data_path, '/Magnetization/error', start, count);
    meanEnergy(:, i) = h5read(data_path, '/Energy/mean', start, count);
    errorEnergy(:, i) = h5read(data_path, '/Energy/error', start, count);
    mean_fee(:, i) = h5read(data_path, '/HeatCapacity/mean', start, count);
    error_fee(:, i) = h5read(data_path, '/HeatCapacity/error', start, count);
    mean_fmm(:, i) = h5read(data_path, '/Susceptibility/mean', start, count);
    error_fmm(:, i) = h5read(data_path, '/Susceptibility/error', start, count);
    mean_fme(:, i) = h5read(data_path, '/mDer/mean', start, count);
    error_fme(:, i) = h5read(data_path, '/mDer/error', start, count);
    meanBinder(:, i) = h5read(data_path, '/Binder/mean', start, count);
    errorBinder(:, i) = h5read(data_path, '/Binder/error', start, count);
    % mean_Qt(:, i) = h5read(data_path, '/Qt/mean', start, count);
    % error_Qt(:, i) = h5read(data_path, '/Qt/error', start, count);
    Temperature(i) = h5read(data_path, '/Temperature');
    SystemLength(i) = h5read(data_path, '/SystemLength');
    SystemSize(i) = h5read(data_path, '/SystemSize');
    for j = 1:Nbin
        mtime1Bins(:, i, j) = h5read(data_path, ['/mtime1Bins/', num2str(j-1)], start, count);
        mtime2Bins(:, i, j) = h5read(data_path, ['/mtime2Bins/', num2str(j-1)], start, count);
        mtime3Bins(:, i, j) = h5read(data_path, ['/mtime3Bins/', num2str(j-1)], start, count);
        mtime4Bins(:, i, j) = h5read(data_path, ['/mtime4Bins/', num2str(j-1)], start, count);
        Etime1Bins(:, i, j) = h5read(data_path, ['/Etime1Bins/', num2str(j-1)], start, count);
        Etime2Bins(:, i, j) = h5read(data_path, ['/Etime2Bins/', num2str(j-1)], start, count);
        Etime3Bins(:, i, j) = h5read(data_path, ['/Etime3Bins/', num2str(j-1)], start, count);
        Etime4Bins(:, i, j) = h5read(data_path, ['/Etime4Bins/', num2str(j-1)], start, count);
        m1E1timeBins(:, i, j) = h5read(data_path, ['/m1E1timeBins/', num2str(j-1)], start, count);
        
    end
end
disp(length(meanMagnetization))


%%---------------------------------------------------------------------------------
%% plot m
fig = figure;
fig.Position  = [0 0 450 300];
clf('reset')

t_max = length(meanMagnetization);
t_range = 1:t_max;
m__rel_error = errorMagnetization ./ meanMagnetization;

for i = 1:length(Temperature)
    % errorbar(t_range, meanMagnetization(t_range,i), errorMagnetization(t_range,i), ...
    %     'Marker', '.', 'MarkerSize', 5);
     % plot(t_range, errorMagnetization(t_range,i), 'Marker', '.', 'MarkerSize', 7);
    plot(t_range, meanMagnetization(t_range,i), 'Marker', '.', 'MarkerSize', 12);
     % plot(t_range, m__rel_error(t_range,i), 'Marker', '.', 'MarkerSize', 7);

    % plot(log(t_range), log(meanMagnetization(t_range,i)), 'Marker', '.', 'MarkerSize', 7);

    hold on
end

% t_range = 1:t_max;
% for i = 1:length(Temperature)/2
% % for i = 2
%     % errorbar(t_range, meanMagnetization(t_range,i), errorMagnetization(t_range,i), ...
%     %     'Marker', '.', 'MarkerSize', 5);
%     % plot(2*4*t_range, errorMagnetization(t_range,i), 'Marker', '.', 'MarkerSize', 7);
%     plot(2*7.5*t_range, meanMagnetization(t_range,i), 'Marker', '.', 'MarkerSize', 7);
%     % plot(log(t_range), log(meanMagnetization(t_range,i)), 'Marker', '.', 'MarkerSize', 7);
% 
%     hold on
% end
% 
% 
% for i = length(Temperature)/2+1:length(Temperature)
% % for i = 2 + 3
%     % errorbar(t_range, meanMagnetization(t_range,i), errorMagnetization(t_range,i), ...
%     %     'Marker', '.', 'MarkerSize', 5);
%     % plot(t_range, errorMagnetization(t_range,i), 'Marker', '.', 'MarkerSize', 7);
%     plot(t_range, meanMagnetization(t_range,i), 'Marker', '.', 'MarkerSize', 7);
%     % plot(log(t_range), log(meanMagnetization(t_range,i)), 'Marker', '.', 'MarkerSize', 7);
% 
%     hold on
% end

% xlim([1 1e3])
m_error_mean = mean(errorMagnetization);
m__rel_error_mean = mean(errorMagnetization ./ meanMagnetization)

xlabel '$t \; [\mathrm{MCS}]$'
ylabel '$|\mathbf{m}|$'
% ylabel '$\mathrm{error}(|\mathbf{m}|)$'

% legendStrings = "$L$ = " + string(SystemLength) + ", $T$ = " + string(Temperature);
% legendStrings = "$T$ = " + string(Temperature);
legendStrings = "$L$ = " + string(SystemLength);

legend(legendStrings, 'Location','northeast', ...
    'Interpreter','latex')
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
% saveas(gcf, [save_address '/', 'm.png']);
% saveas(gcf, [save_address '/', 'm_scaled.png']);
exportgraphics(gcf,[save_address '/', 'm_scaled.png'],'Resolution',1200)
% saveas(gcf, [save_address '/', 'merror.png']);
% saveas(gcf, [save_address '/', 'merror_scaled.png']);
%%---------------------------------------------------------------------------------
%% plot E
clf('reset')
t_max = length(meanMagnetization);
t_range = 1:max_count;

for i = 1:length(Temperature)
    % errorbar(t_range, meanMagnetization(t_range,i), errorMagnetization(t_range,i), ...
    %     'Marker', '.', 'MarkerSize', 5);
    plot(t_range, meanEnergy(t_range,i), 'Marker', '.', 'MarkerSize', 7);

    hold on
end

% t_range = 1:t_max;
% for i = 1:length(Temperature)/2
% % for i = 2
%     % errorbar(t_range, meanMagnetization(t_range,i), errorMagnetization(t_range,i), ...
%     %     'Marker', '.', 'MarkerSize', 5);
%     % plot(2*4*t_range, errorEnergy(t_range,i), 'Marker', '.', 'MarkerSize', 7);
%     plot(2*7.5*t_range, meanEnergy(t_range,i), 'Marker', '.', 'MarkerSize', 7);
%     % plot(log(t_range), log(meanMagnetization(t_range,i)), 'Marker', '.', 'MarkerSize', 7);
% 
%     hold on
% end
% 
% 
% for i = length(Temperature)/2+1:length(Temperature)
% % for i = 2 + 3
%     % errorbar(t_range, meanMagnetization(t_range,i), errorMagnetization(t_range,i), ...
%     %     'Marker', '.', 'MarkerSize', 5);
%     % plot(t_range, errorEnergy(t_range,i), 'Marker', '.', 'MarkerSize', 7);
%     plot(t_range, meanEnergy(t_range,i), 'Marker', '.', 'MarkerSize', 7);
%     % plot(log(t_range), log(meanMagnetization(t_range,i)), 'Marker', '.', 'MarkerSize', 7);
% 
%     hold on
% end

E_error_mean = mean(errorEnergy)

xlabel '$t$'
ylabel '$E$'
legendStrings = "$L$ = " + string(SystemLength) + ", $T$ = " + string(Temperature);
% legendStrings = "$T$ = " + string(Temperature);
legend(legendStrings, 'Location','southeast', ...
    'Interpreter','latex')
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
% saveas(gcf, [save_address '/', 'E.png']);
saveas(gcf, [save_address '/', 'E_scaled.png']);
%%---------------------------------------------------------------------------------
%% Delta M with larger size

index_of_reference = 2;
errorDeltaM = [];
TotalSizes = length(Temperature);
SystemRange = [1:index_of_reference-1 index_of_reference+1:TotalSizes];
for i = 1:size(meanMagnetization, 2)
    data1 = meanMagnetization(:, i);
    data2 = meanMagnetization(:, index_of_reference);
    
    error1 = errorMagnetization(:, i);
    error2 = errorMagnetization(:, index_of_reference);
    % Calculate the absolute differences between data points
    diff = abs(data1 - data2);
    
    % Calculate the sum of error bars
    errorSum = error1 + error2;
    
    % Find where the difference is less than or equal to the sum of error bars
    equalIndices = diff <= errorSum;
    
    % Initialize the result vector with the minimum distances
    errorDeltaM(:, i) = (diff - errorSum) ./ data2;
    
    % Set the indices where data points are equal within error bars to 0
    errorDeltaM(equalIndices, i) = 0;
end

clf('reset')
t_max = length(meanMagnetization);
t_range = 1:1:t_max;
for i =  SystemRange
    plot(t_range, errorDeltaM(t_range,i), 'Marker', '.', 'MarkerSize', 5);
    hold on
end

xlabel '$t$'
ylabel '$\Delta m / m$'
legendStrings = "$L$ = " + string(SystemLength(SystemRange)) + ", $T$ = " + string(Temperature(SystemRange));
% legendStrings = "$T$ = " + string(Temperature);
legend(legendStrings, 'Location','best', ...
    'Interpreter','latex','fontsize',12)
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'DeltaM.png']);

%%---------------------------------------------------------------------------------
%% Delta E with larger size

index_of_reference = 1;
errorDeltaM = [];
for i = 1:size(meanEnergy, 2)
    data1 = meanEnergy(:, i);
    data2 = meanEnergy(:, index_of_reference);
    
    error1 = errorEnergy(:, i);
    error2 = errorEnergy(:, index_of_reference);
    % Calculate the absolute differences between data points
    diff = abs(data1 - data2);
    
    % Calculate the sum of error bars
    errorSum = error1 + error2;
    
    % Find where the difference is less than or equal to the sum of error bars
    equalIndices = diff <= errorSum;
    
    % Initialize the result vector with the minimum distances
    errorDeltaM(:, i) = (diff - errorSum) ./ abs(data2);
    
    % Set the indices where data points are equal within error bars to 0
    errorDeltaM(equalIndices, i) = 0;
end

clf('reset')
t_max = length(meanEnergy);
t_range = 1:1:1000;
for i = 1:length(Temperature)
    plot(t_range, errorDeltaM(t_range,i), 'Marker', '.', 'MarkerSize', 5);
    hold on
end

xlabel '$t$'
ylabel '$\Delta E / |E|$'
legendStrings = "$L$ = " + string(SystemLength) + ", $T$ = " + string(Temperature);
% legendStrings = "$T$ = " + string(Temperature);
legend(legendStrings, 'Location','northeast', ...
    'Interpreter','latex','fontsize',12)
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'DeltaE.png']);
%%---------------------------------------------------------------------------------
%% plot fmm
clf('reset')
t_max = length(meanMagnetization);
t_range = 1:t_max;

for i = 1:length(Temperature)
    errorbar(t_range, mean_fmm(t_range,i), error_fmm(t_range,i), ...
        'Marker', '.', 'MarkerSize', 5);
    hold on
end


xlabel '$t$'
ylabel '$f_{mm}$'
legendStrings = "$L$ = " + string(SystemLength) + ", $T$ = " + string(Temperature);
% legendStrings = "$T$ = " + string(Temperature);
legend(legendStrings, 'Location','southeast', ...
    'Interpreter','latex','fontsize',12)
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'fmm.png']);

%%---------------------------------------------------------------------------------
%% plot fmm * L2
clf('reset')
t_max = length(meanMagnetization);
t_range = 1:t_max;
A = [];
for i = 1:length(Temperature)
    A(:, i) = mean_fmm(t_range,i) .* double(SystemLength(i)) ^ 2;
    Aerr(:, i) = error_fmm(t_range,i) .* double(SystemLength(i)) ^ 2;
    errorbar(t_range, A(:, i), Aerr(:, i), 'Marker', '.', 'MarkerSize', 10);
    hold on
end

errMean = mean(Aerr)

xlabel '$t$'
ylabel '$f_{mm} L^2$'
legendStrings = "$L$ = " + string(SystemLength) + ", $T$ = " + string(Temperature);
% legendStrings = "$T$ = " + string(Temperature);
legend(legendStrings, 'Location','southeast', ...
    'Interpreter','latex','fontsize',12)
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'fmmL2.png']);

%%---------------------------------------------------------------------------------
%% plot fme
clf('reset')
t_max = length(meanMagnetization);
t_range = 1:t_max;

for i = 1:length(Temperature)
    errorbar(t_range, mean_fme(t_range,i), error_fme(t_range,i), ...
        'Marker', '.', 'MarkerSize', 5);
    % plot(t_range, mean_fme(t_range,i), 'Marker', '.', 'MarkerSize', 5);
    hold on
end

xlabel '$t$'
ylabel '$f_{me}$'
legendStrings = "$L$ = " + string(SystemLength) + ", $T$ = " + string(Temperature);
% legendStrings = "$T$ = " + string(Temperature);
legend(legendStrings, 'Location','southeast', ...
    'Interpreter','latex','fontsize',12)
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'fme.png']);


%%---------------------------------------------------------------------------------
%% plot fme * L2
clf('reset')
t_max = length(meanMagnetization);
t_range = 1:t_max;
A = [];
for i = 1:length(Temperature)
    A(:, i) = mean_fme(t_range,i) .* double(SystemLength(i)) ^ 2;
    Aerr(:, i) = error_fme(t_range,i) .* double(SystemLength(i)) ^ 2;
    errorbar(t_range, A(:, i), Aerr(:, i), 'Marker', '.', 'MarkerSize', 10);
    hold on
end

AerrMean = mean(Aerr)

xlabel '$t$'
ylabel '$f_{me} L^2$'
legendStrings = "$L$ = " + string(SystemLength) + ", $T$ = " + string(Temperature);
% legendStrings = "$T$ = " + string(Temperature);
legend(legendStrings, 'Location','southeast', ...
    'Interpreter','latex','fontsize',12)
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'fmeL2.png']);

%%---------------------------------------------------------------------------------
%% plot fee
clf('reset')
t_max = length(meanMagnetization);
t_range = 1:t_max;

for i = 1:length(Temperature)
    errorbar(t_range, mean_fee(t_range,i), error_fee(t_range,i), ...
        'Marker', '.', 'MarkerSize', 5);
    hold on
end

xlabel '$t$'
ylabel '$f_{ee}$'
legendStrings = "$L$ = " + string(SystemLength) + ", $T$ = " + string(Temperature);
% legendStrings = "$T$ = " + string(Temperature);
legend(legendStrings, 'Location','southeast', ...
    'Interpreter','latex','fontsize',12)
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'fee.png']);

%%---------------------------------------------------------------------------------
%% plot fee * L2
clf('reset')
t_max = length(meanMagnetization);
t_range = 1:t_max;
A = [];
for i = 1:length(Temperature)
    A(:, i) = mean_fee(t_range,i) .* double(SystemLength(i)) ^ 2;
    Aerr(:, i) = error_fee(t_range,i) .* double(SystemLength(i)) ^ 2;
    errorbar(t_range, A(:, i), Aerr(:, i), 'Marker', '.', 'MarkerSize', 10);
    hold on
end

AerrMean = mean(Aerr)

xlabel '$t$'
ylabel '$f_{ee} L^2$'
legendStrings = "$L$ = " + string(SystemLength) + ", $T$ = " + string(Temperature);
% legendStrings = "$T$ = " + string(Temperature);
legend(legendStrings, 'Location','southeast', ...
    'Interpreter','latex','fontsize',12)
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'feeL2.png']);

%%---------------------------------------------------------------------------------
%% plot Qt
clf('reset')
t_max = length(meanMagnetization);
t_range = 1:t_max;

for i = 1:length(Temperature)
    errorbar(t_range, mean_Qt(t_range,i), error_Qt(t_range,i), ...
        'Marker', '.', 'MarkerSize', 5);
    hold on
end

xlabel '$t$'
ylabel '$Q(t)$'
legendStrings = "$L$ = " + string(SystemLength) + ", $T$ = " + string(Temperature);
% legendStrings = "$T$ = " + string(Temperature);
legend(legendStrings, 'Location','southeast', ...
    'Interpreter','latex','fontsize',12)
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'fme.png']);

%%---------------------------------------------------------------------------------
%% plot Rushbrooke
clf('reset')
t_max = length(meanMagnetization);
t_range = 1:t_max;

for i = 1:length(Temperature)
    R = (mean_fmm(t_range,i) .* mean_fee(t_range,i)) ./ (mean_fme(t_range,i) .^ 2);
    plot(t_range, R,'Marker', '.', 'MarkerSize', 10);
    hold on
end

xlabel '$t$'
ylabel 'Rushbrooke function'
legendStrings = "$L$ = " + string(SystemLength) + ", $T$ = " + string(Temperature);
% legendStrings = "$T$ = " + string(Temperature);
legend(legendStrings, 'Location','southeast', ...
    'Interpreter','latex','fontsize',12)
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'Rushbrooke.png']);


%%---------------------------------------------------------------------------------
%% plot dTaulnM

% use Finite Difference Formulae for Unequal Sub-Intervals Using Lagrange's Interpolation Formula
clf('reset')
t_max = length(meanMagnetization);
t_range = 1:t_max;

i_before = 2;
i_calc = 3;
i_after = 4;
h1 = Temperature(i_calc) - Temperature(i_before);
h2 = Temperature(i_after) - Temperature(i_calc);

% dTaulnM = - (log(meanMagnetization(:, i_after)) - log(meanMagnetization(:, i_before)))...
%             / (h1 + h2);

f0 = log(meanMagnetization(:, i_before));
f1 = log(meanMagnetization(:, i_calc));
f2 = log(meanMagnetization(:, i_after));
dTaulnM = - (- h2/(h1*(h1+h2)) * f0 - (h1-h2)/(h1*h2) * f1 + h1/(h2*(h1+h2)) * f2);

% h = Temperature(4) - Temperature(1);
% dTaulnM = - (log(meanMagnetization(:, 4)) - log(meanMagnetization(:, 1)))...
%             / (2 * h);


plot(t_range, dTaulnM(t_range), 'Marker', '.', 'MarkerSize', 8);
hold on

t_range = 70:t_max;
logt = log(t_range);
logt = logt';
logm = log(dTaulnM(t_range));

% Fit the linear model
A = logm;
mdl = fitlm(logt, A , 'linear');

% Get the R-squared value
r_squared = mdl.Rsquared.Ordinary;
R2 = r_squared;
slope = mdl.Coefficients.Estimate(2);
error = mdl.Coefficients.SE(2);

disp(['R2 = ' num2str(R2)]);
disp(['slope = ' fmtMeanUnc(slope, error)]);
% slope_max = slope(1)
axis padded




xlabel '$t$'
ylabel '$-\partial_\tau \ln(M)$'
legendStrings = "$L$ = " + string(SystemLength(i_calc)) + ", $T$ = " + string(Temperature(i_calc));
% legendStrings = "$T$ = " + string(Temperature);
legend(legendStrings, 'Location','southeast', ...
    'Interpreter','latex','fontsize',12)
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'dTaulnM.png']);

%%---------------------------------------------------------------------------------
%% find 1/(nu*z) using dTaulnM with errorbar
t_max = max_count;
t_range = 30:t_max;
logt = log(t_range);
logt = logt';

i_before = 1;
i_calc = 2;
i_after = 3;
h1 = Temperature(i_calc) - Temperature(i_before);
h2 = Temperature(i_after) - Temperature(i_calc);

invNuZ_Bin = [];

for Bi = 1:5
    m1 = squeeze(mtime1Bins(t_range, i_calc, Bi));
    m2 = squeeze(mtime2Bins(t_range, i_calc, Bi));
    m3 = squeeze(mtime3Bins(t_range, i_calc, Bi));
    m4 = squeeze(mtime4Bins(t_range, i_calc, Bi));
    E1 = squeeze(Etime1Bins(t_range, i_calc, Bi));
    E2 = squeeze(Etime2Bins(t_range, i_calc, Bi));
    E3 = squeeze(Etime3Bins(t_range, i_calc, Bi));
    E4 = squeeze(Etime4Bins(t_range, i_calc, Bi));
    m1E1 = squeeze(m1E1timeBins(t_range, i_calc, Bi));


    f0 = squeeze(mtime1Bins(t_range, i_before, Bi));
    f1 = squeeze(mtime1Bins(t_range, i_calc, Bi));
    f2 = squeeze(mtime1Bins(t_range, i_after, Bi));
    dTaulnM = - (- h2/(h1*(h1+h2)) * f0 - (h1-h2)/(h1*h2) * f1 + h1/(h2*(h1+h2)) * f2);

    p = polyfit(logt, log(dTaulnM), 1);
    invNuZ_Bin(Bi) = p(1);
end

disp(['1/(nu * z) = ' computeAndFormatMean(invNuZ_Bin), newline]);
%%---------------------------------------------------------------------------------
%% plot R2
clf('reset')
t_max = max_count;
t_range = 100:t_max;
logt = log(t_range);
logt = logt';
logm = log(meanMagnetization(t_range,:));
% logm = log(mean_fmm(t_range,:));
% logm = log(mean_fme(t_range,:));
% logm = log(mean_fee(t_range,:));
% logm = log(mean_Qt(t_range,:));
% logt_smooth = smooth(log(t_range), 100);

R2Vec = zeros(size(Temperature));
slope = zeros(size(Temperature));
slopeErr = zeros(size(Temperature));
for i = 1:length(Temperature)
    % Fit the linear model
    A = logm(:,i);
    mdl = fitlm(logt, A , 'linear');
    
    % Get the R-squared value
    r_squared = mdl.Rsquared.Ordinary;
    R2Vec(i) = r_squared;
    slope(i) = mdl.Coefficients.Estimate(2);
    slopeErr(i) = mdl.Coefficients.SE(2);
end

[~, idx] = max(R2Vec);
% errorbar(Temperature, slope, slopeErr, '.','MarkerSize', 15 ); xlabel '$T$'; ylabel '$slope$'
% errorbar(SystemLength, slope, slopeErr, '.','MarkerSize', 15 ); xlabel '$T$'; ylabel '$slope$'
plot(Temperature, R2Vec, 'Marker', '.','MarkerSize', 15 ); xlabel '$T$'; ylabel '$R^2$'
% plot(SystemLength, R2Vec, '.','MarkerSize', 15 ); xlabel '$L$'; ylabel '$R^2$'
% 
% x = Temperature;
% y = R2Vec;
% degree = 2; % Change this to the degree of the polynomial you want
% coefficients = polyfit(x, y, degree);
% 
% x_fit = linspace(min(x), max(x), 500); % Adjust the range accordingly
% y_fit = polyval(coefficients, x_fit);
% 
% plot(Temperature, R2Vec, '.', x_fit, y_fit, '-', ...
%     'MarkerSize', 15, 'LineWidth', 1)
% legend('Original Data', 'Parabolic Fit')
% idx
% slope
% slopeErr
% slope_max = slope(1)
axis padded

% yl = ylim;
% xl = xlim;
% yrange = (yl(2)- yl(1));
% xrange = (xl(2)- xl(1))
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));



set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
% cmap = linspecer(length(Temperature));
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
% saveas(gcf, [save_address '/', 'm.png']);
% saveas(gcf, [save_address '/', 'slope.png']);
saveas(gcf, [save_address '/', 'R2.png']);

%%---------------------------------------------------------------------------------
%% *** find critical exponents with errorbar
t_max = max_count;
t_range = 50:500;
logt = log(t_range);
logt = logt';
i_calc = 3;
Nboot = 100;
zBin = [];
nuBin = [];
betaBin = [];
alphaBin = [];

lambda_mBin = [];
lambda_mmBin = [];
lambda_meBin = [];
lambda_eeBin = [];

lambda_m_max_fit_error = 0;
lambda_mm_max_fit_error = 0;
lambda_me_max_fit_error = 0;
lambda_ee_max_fit_error = 0;


for Bi = 1:Nboot
    BootBins = randi(Nbin, Nbin, 1);
    BootBins = BootBins';
    m1 = BootAverage(mtime1Bins, BootBins);
    m2 = BootAverage(mtime2Bins, BootBins);
    m3 = BootAverage(mtime3Bins, BootBins);
    m4 = BootAverage(mtime4Bins, BootBins);
    E1 = BootAverage(Etime1Bins, BootBins);
    E2 = BootAverage(Etime2Bins, BootBins);
    E3 = BootAverage(Etime3Bins, BootBins);
    E4 = BootAverage(Etime4Bins, BootBins);
    m1E1 = BootAverage(m1E1timeBins, BootBins);

    m1 = squeeze(m1(t_range, i_calc));
    m2 = squeeze(m2(t_range, i_calc));
    m3 = squeeze(m3(t_range, i_calc));
    m4 = squeeze(m4(t_range, i_calc));
    E1 = squeeze(E1(t_range, i_calc));
    E2 = squeeze(E2(t_range, i_calc));
    E3 = squeeze(E3(t_range, i_calc));
    E4 = squeeze(E4(t_range, i_calc));
    m1E1 = squeeze(m1E1(t_range, i_calc));

    fmm = m2 ./ (m1 .^2) - 1;
    fee = E2 ./ (E1 .^2) - 1;
    fme = m1E1 ./ (m1 .* E1) - 1;

    % p = polyfit(logt, log(m1), 1);
    % lambda_m = p(1);
    % 
    % p = polyfit(logt, log(fmmBin), 1);
    % lambda_mm = p(1);
    % 
    % p = polyfit(logt, log(fmeBin), 1);
    % lambda_me = p(1);
    % 
    % p = polyfit(logt, log(feeBin), 1);
    % lambda_ee = p(1);

    mdl = fitlm(logt, log(m1) , 'linear');
    lambda_m = mdl.Coefficients.Estimate(2);
    lambda_m_max_fit_error = max(lambda_m_max_fit_error, mdl.Coefficients.SE(2));

    mdl = fitlm(logt, log(fmm) , 'linear');
    lambda_mm = mdl.Coefficients.Estimate(2);
    lambda_mm_max_fit_error = max(lambda_mm_max_fit_error, mdl.Coefficients.SE(2));
    
    mdl = fitlm(logt, log(fme) , 'linear');
    lambda_me = mdl.Coefficients.Estimate(2);
    lambda_me_max_fit_error = max(lambda_me_max_fit_error, mdl.Coefficients.SE(2));

    mdl = fitlm(logt, log(fee) , 'linear');
    lambda_ee = mdl.Coefficients.Estimate(2);
    lambda_ee_max_fit_error = max(lambda_ee_max_fit_error, mdl.Coefficients.SE(2));

    zBoot(Bi) = 2 / lambda_mm;
    nuBoot(Bi) = lambda_mm / (2 * lambda_me);
    betaBoot(Bi) = lambda_m / lambda_me;
    alphaBoot(Bi) = lambda_ee / lambda_me;
    

    lambda_mBoot(Bi) = lambda_m;
    lambda_mmBoot(Bi) = lambda_mm;
    lambda_meBoot(Bi) = lambda_me;
    lambda_eeBoot(Bi) = lambda_ee;
end
z_mean = mean(zBoot);
nu_mean = mean(nuBoot);
beta_mean = mean(betaBoot);
z_error = std(zBoot);
nu_error = std(nuBoot);
beta_error = std(betaBoot);

disp(['z = ' fmtMeanUnc(z_mean, z_error)]);
disp(['nu = ' fmtMeanUnc(nu_mean, nu_error)]);
disp(['beta = ' fmtMeanUnc(beta_mean, beta_error)]);
% disp(['alpha = ' computeAndFormatMean(alphaBin)]);
disp(['beta/nu = ' num2str(beta_mean/nu_mean)]);
% disp(['eta = ' computeAndFormatMean(etaBin), newline]);


% disp(['lambda_m = ' computeAndFormatMean(lambda_mBin)]);
% disp(['lambda_mm = ' computeAndFormatMean(lambda_mmBin)]);
% disp(['lambda_me = ' computeAndFormatMean(lambda_meBin)]);
% disp(['lambda_ee = ' computeAndFormatMean(lambda_eeBin), newline]);

% for Bi = 1:5
%     logt = log(t_range);
%     logt = logt';
%     m1 = squeeze(mtime1Bins(t_range, i_calc, Bi));
%     m2 = squeeze(mtime2Bins(t_range, i_calc, Bi));
%     m3 = squeeze(mtime3Bins(t_range, i_calc, Bi));
%     m4 = squeeze(mtime4Bins(t_range, i_calc, Bi));
%     E1 = squeeze(Etime1Bins(t_range, i_calc, Bi));
%     E2 = squeeze(Etime2Bins(t_range, i_calc, Bi));
%     E3 = squeeze(Etime3Bins(t_range, i_calc, Bi));
%     E4 = squeeze(Etime4Bins(t_range, i_calc, Bi));
%     m1E1 = squeeze(m1E1timeBins(t_range, i_calc, Bi));
% 
%     fmmBin = m2 ./ (m1 .^2) - 1;
%     feeBin = E2 ./ (E1 .^2) - 1;
%     fmeBin = m1E1 ./ (m1 .* E1) - 1;
% 
%     % p = polyfit(logt, log(m1), 1);
%     % lambda_m = p(1);
%     % 
%     % p = polyfit(logt, log(fmmBin), 1);
%     % lambda_mm = p(1);
%     % 
%     % p = polyfit(logt, log(fmeBin), 1);
%     % lambda_me = p(1);
%     % 
%     % p = polyfit(logt, log(feeBin), 1);
%     % lambda_ee = p(1);
% 
%     mdl = fitlm(logt, log(m1) , 'linear');
%     lambda_m = mdl.Coefficients.Estimate(2);
%     lambda_m_max_fit_error = max(lambda_m_max_fit_error, mdl.Coefficients.SE(2));
% 
%     mdl = fitlm(logt, log(fmmBin) , 'linear');
%     lambda_mm = mdl.Coefficients.Estimate(2);
%     lambda_mm_max_fit_error = max(lambda_mm_max_fit_error, mdl.Coefficients.SE(2));
% 
%     mdl = fitlm(logt, log(fmeBin) , 'linear');
%     lambda_me = mdl.Coefficients.Estimate(2);
%     lambda_me_max_fit_error = max(lambda_me_max_fit_error, mdl.Coefficients.SE(2));
% 
%     mdl = fitlm(logt, log(feeBin) , 'linear');
%     lambda_ee = mdl.Coefficients.Estimate(2);
%     lambda_ee_max_fit_error = max(lambda_ee_max_fit_error, mdl.Coefficients.SE(2));
% 
%     zBin(Bi) = 2 / lambda_mm;
%     nuBin(Bi) = lambda_mm / (2 * lambda_me);
%     betaBin(Bi) = lambda_m / lambda_me;
%     alphaBin(Bi) = lambda_ee / lambda_me;
% 
% 
%     lambda_mBin(Bi) = lambda_m;
%     lambda_mmBin(Bi) = lambda_mm;
%     lambda_meBin(Bi) = lambda_me;
%     lambda_eeBin(Bi) = lambda_ee;
% end
% z_mean = mean(zBin);
% nu_mean = mean(nuBin);
% beta_mean = mean(betaBin);
% z_error = std(zBin, 0) / sqrt(length(zBin));
% nu_error = std(nuBin, 0) / sqrt(length(nuBin));
% beta_error = std(betaBin, 0) / sqrt(length(betaBin));
% 
% 
% 
% disp(['z = ' fmtMeanUnc(z_mean, z_error)]);
% disp(['nu = ' fmtMeanUnc(nu_mean, nu_error)]);
% disp(['beta = ' fmtMeanUnc(beta_mean, beta_error)]);
% disp(['alpha = ' computeAndFormatMean(alphaBin)]);
% disp(['beta/nu = ' num2str(beta_mean/nu_mean)]);
% % disp(['eta = ' computeAndFormatMean(etaBin), newline]);
% 
% disp(['lambda_m = ' computeAndFormatMean(lambda_mBin)]);
% disp(['lambda_mm = ' computeAndFormatMean(lambda_mmBin)]);
% disp(['lambda_me = ' computeAndFormatMean(lambda_meBin)]);
% disp(['lambda_ee = ' computeAndFormatMean(lambda_eeBin), newline]);






% 
% lambda_m_max_fit_error
% lambda_mm_max_fit_error
% lambda_me_max_fit_error
% lambda_ee_max_fit_error
%%---------------------------------------------------------------------------------
%% find critical exponents
t_max = max_count;
t_range = 100:1000;
logt = log(t_range);
logt = logt';
i_calc = 5;

logt = log(t_range);
logt = logt';

m1 = meanMagnetization(t_range,i_calc);
fmm = mean_fmm(t_range,i_calc);
fme = mean_fme(t_range,i_calc);
fee = mean_fee(t_range,i_calc);

mdl = fitlm(logt, log(m1) , 'linear');
lambda_m = mdl.Coefficients.Estimate(2);

mdl = fitlm(logt, log(fmm) , 'linear');
lambda_mm = mdl.Coefficients.Estimate(2);

mdl = fitlm(logt, log(fme) , 'linear');
lambda_me = mdl.Coefficients.Estimate(2);

mdl = fitlm(logt, log(fee) , 'linear');
lambda_ee = mdl.Coefficients.Estimate(2);


z = 2 / lambda_mm;
nu = lambda_mm / (2 * lambda_me);
beta = lambda_m / lambda_me;
alpha = lambda_ee / lambda_me;

disp(['z = ' num2str(z)]);
disp(['nu = ' num2str(nu)]);
disp(['beta = ' num2str(beta)]);
% disp(['alpha = ' num2str(alpha)]);
disp(['beta/nu = ' num2str(beta/nu), newline]);
% disp(['eta = ' computeAndFormatMean(etaBin), newline]);

% disp(['lambda_m = ' computeAndFormatMean(lambda_mBin)]);
% disp(['lambda_mm = ' computeAndFormatMean(lambda_mmBin)]);
% disp(['lambda_me = ' computeAndFormatMean(lambda_meBin)]);
% disp(['lambda_ee = ' computeAndFormatMean(lambda_eeBin), newline]);
%%---------------------------------------------------------------------------------
%% plot trend
t_max = 1000;
t_range = 1:t_max;
logt = log(t_range);
logt = logt';
logm = log(meanMagnetization(t_range,:));
logt_smooth = smooth(log(t_range), 100);
divt = 1 ./ t_range;
NTemp = length(Temperature);
clf('reset')
grad = [];
average_grad= [];
trend = [];
AverageStep = 100;
for i = 1:NTemp
    logm_smooth = smooth(log(meanMagnetization(t_range,i)), 100);
    grad(:,i) = - gradient(logm_smooth, logt_smooth);
    ones_vector = ones(1, AverageStep);
    conv_grad = conv(grad(:,i), ones_vector, 'valid');
    average_grad(:,i) = conv_grad / AverageStep;
    conv_time = conv(t_range, ones_vector, 'valid');
    average_time = conv_time / AverageStep;
    
    data = logm(:, i);
    n = length(data);
    % derivative(:, i) = zeros(1, n);
    for j = AverageStep+1:n
        p = polyfit(logt(j-AverageStep:j), - data(j-AverageStep:j), 1);
        trend(j, i) = p(1);
        % trend(j, i) = - (data(j-AverageStep) - data(j)) / (logt(j-AverageStep) - logt(j));
    end
end

t_range_final = (AverageStep+71):1000;
for i = 1:NTemp
    plot(divt(t_range_final), trend(t_range_final,i) ,'Marker', '.','MarkerSize', 15);
    hold on
end

axis normal

% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
xlabel '$t^{-1}$'
ylabel 'Trend in the slope of $\log(|\mathbf{m}|)$'

legendStrings = "$L$ = " + string(SystemLength) + ", $T$ = " + string(Temperature);
% legendStrings = "$T$ = " + string(Temperature);
legend(legendStrings, 'Location','southeast','Interpreter','latex','fontsize',12)
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(NTemp);
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'trend.png']);

%%---------------------------------------------------------------------------------
%% points above below proportion


t_max = max_count;
t_range = 70:t_max;
logt = log(t_range);
logt = logt';
NTemp = length(Temperature);
i = 1;
logm_smooth = smooth(log(meanMagnetization(t_range,i)), 100);
logm_unsmooth = log(meanMagnetization(t_range,i));

x = logt;
y = logm_unsmooth;
% Start and end points of the curve
x_start = x(1);
y_start = y(1);
x_end = x(end);
y_end = y(end);

% Equation of the line from the start to the end of the curve
a = (y_end - y_start) / (x_end - x_start); % slope of the line
b = y_start - a * x_start; % y-intercept of the line

% Line's Y values across the curve range
line_y = a * x + b;

% Logical indices of points above and below the line
points_above = y > line_y;
points_below = y < line_y;

% Proportion of points
proportion_above = sum(points_above) / length(y);
proportion_below = sum(points_below) / length(y);

% Output results
fprintf('Proportion above the line: %f\n', proportion_above);
fprintf('Proportion below the line: %f\n', proportion_below);

% Optional: plotting the curve and the line
plot(x, y, 'b', x, line_y, 'r--');
legend('Curve', 'Line from start to end');
xlabel('X');
ylabel('Y');
title('Curve and Line');

%%---------------------------------------------------------------------------------
%% points above proportion multiple
t_max = max_count;
t_range = 100:1000;
logt = log(t_range);
logt = logt';
NTemp = length(Temperature);

% cmap = linspecer(NTemp);
cmap = haxby(NTemp+1);
mi = 1;

% Make sure that x_cell and y_cell are the same length
numCurves = length(Temperature);

% Initialize arrays to store proportions for each curve
proportion_above = zeros(1, numCurves);
proportion_below = zeros(1, numCurves);

% Create figure outside of loop to hold all plots
% figure;
hold on; % Allows multiple plots in the same figure

% Process each curve
mi = 1;
for i = 1:numCurves
    % Extract the current curve
    logm_smooth = smooth(log(meanMagnetization(t_range,i)), 100);
    logm_unsmooth = log(meanMagnetization(t_range,i));
    x = logt;
    y = logm_unsmooth;
    
    % Start and end points of the curve
    x_start = x(1);
    y_start = y(1);
    x_end = x(end);
    y_end = y(end);
    
    % Equation of the line
    a = (y_end - y_start) / (x_end - x_start); % slope of the line
    b = y_start - a * x_start; % y-intercept
    
    % Calculate the line values
    line_y = a * x + b;
    
    % Determine points above and below the line
    points_above = y > line_y;
    points_below = y < line_y;
    
    % Calculate proportions
    proportion_above(i) = sum(points_above) / length(y);
    proportion_below(i) = sum(points_below) / length(y);
    
    % Plot the curve and its line on the same graph
    h = plot(x, y,...
     'Color',cmap(mi,:), 'Marker', '.','MarkerSize', 1, 'LineWidth', 2, ...
     'LineStyle', '-' , 'DisplayName',  ['T = ', num2str(Temperature(i))]);
    plot(x, line_y, '--', 'HandleVisibility', 'off','Color',cmap(mi,:)); 
    mi = mi + 1;
end

% Add graph details
xlabel '$\log(t)$'
ylabel '$$\log(m)$'
legendStrings = "$L$ = " + string(SystemLength) + ", $T$ = " + string(Temperature);
legend(legendStrings, 'Location','best', ...
    'Interpreter','latex')
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

figure;
plot(Temperature, proportion_above, 'Marker', 'o'); 
xlabel '$T$'
ylabel 'Proportion above line'
ylim([0 1])
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
% Output results
for i = 1:numCurves
    % fprintf('%0.3f - above the line: %0.3f\n', Temperature(i), proportion_above(i));
    fprintf('T = %0.4f : %0.3f proportion of points above line\n', Temperature(i), proportion_above(i));
end
fprintf('\n')
saveas(gcf, [save_address '/', 'points above line.png']);

%%---------------------------------------------------------------------------------
%% points above proportion multiple errorbar
fig = figure;
fig.Position  = [0 0 300 300];
t_max = max_count;
t_range = 100:1000;
Nboot = 100;
NTemp = length(Temperature);
pVec = zeros(NTemp, Nboot);
for iboot = 1:Nboot
    BootBins = randi(Nbin, Nbin, 1);
    BootBins = BootBins';
    logt = log(t_range);
    logt = logt';
    for i = 1:NTemp
        mboot = zeros(size(squeeze(mtime1Bins(t_range, i, 1))));
        for Bi = BootBins
            mboot = mboot + squeeze(mtime1Bins(t_range, i, Bi));
        end
        mboot = mboot / Nbin;
        logm= log(mboot);

        x = logt;
        y = logm;

        % Start and end points of the curve
        x_start = x(1);
        y_start = y(1);
        x_end = x(end);
        y_end = y(end);

        % Equation of the line
        a = (y_end - y_start) / (x_end - x_start); % slope of the line
        b = y_start - a * x_start; % y-intercept

        % Calculate the line values
        line_y = a * x + b;

        % Determine points above and below the line
        points_above = y > line_y;
        points_below = y < line_y;

        % Calculate proportions
        proportion_above = sum(points_above) / length(y);
        proportion_below = sum(points_below) / length(y);

        pVec(i, iboot) = proportion_above;
    end
end

pVecMean = mean(pVec,2);
pVecError = std(pVec, 0, 2);

% errorbar(SystemLength, pVecMean, pVecError, '.','MarkerSize', 15 ); 
errorbar(Temperature, pVecMean, pVecError, ...
    'Marker', '.','MarkerSize', 20, 'LineWidth', 0.6, 'Color', 'r' ); 

% line(xlim(), [0,0], 'LineWidth', 0.5, 'Color', 'k');
% grid on;
axis padded
ylim([0 1])
xlabel '$T$'; 
ylabel 'Proportion above line'
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
% cmap = linspecer(length(Temperature));
% ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'points above line errorbar.png']);

%%---------------------------------------------------------------------------------
%% points below proportion multiple errorbar
fig = figure;
fig.Position  = [0 0 450 300];
t_max = max_count;
t_range = 100:1000;
Nboot = 500;
NTemp = length(Temperature);
pVec = zeros(NTemp, Nboot);
for iboot = 1:Nboot
    BootBins = randi(Nbin, Nbin, 1);
    BootBins = BootBins';
    logt = log(t_range);
    logt = logt';
    for i = 1:NTemp
        mboot = zeros(size(squeeze(mtime1Bins(t_range, i, 1))));
        for Bi = BootBins
            mboot = mboot + squeeze(mtime1Bins(t_range, i, Bi));
        end
        mboot = mboot / Nbin;
        logm= log(mboot);

        x = logt;
        y = logm;

        % Start and end points of the curve
        x_start = x(1);
        y_start = y(1);
        x_end = x(end);
        y_end = y(end);

        % Equation of the line
        a = (y_end - y_start) / (x_end - x_start); % slope of the line
        b = y_start - a * x_start; % y-intercept

        % Calculate the line values
        line_y = a * x + b;

        % Determine points above and below the line
        points_above = y > line_y;
        points_below = y < line_y;

        % Calculate proportions
        proportion_above = sum(points_above) / length(y);
        proportion_below = sum(points_below) / length(y);

        pVec(i, iboot) = proportion_below;
    end
end

pVecMean = mean(pVec,2);
pVecError = std(pVec, 0, 2);

colororder({'b','r'})
yyaxis left

% errorbar(SystemLength, pVecMean, pVecError, '.','MarkerSize', 15 ); 
errorbar(Temperature, pVecMean, pVecError, ...
    'Marker', '.','MarkerSize', 20, 'LineWidth', 0.6, 'Color', 'b' ); 

% line(xlim(), [0,0], 'LineWidth', 0.5, 'Color', 'k');
% grid on;
axis padded
ylim([0 1])
xlabel '$T$'; 
ylabel 'Proportion below line'
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
% cmap = linspecer(length(Temperature));
% ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
% saveas(gcf, [save_address '/', 'points below line errorbar.png']);
%%---------------------------------------------------------------------------------
%% qubic curvature
t_max = max_count;
t_range = 100:1000;
logt = log(t_range);
logt = logt';
NTemp = length(Temperature);
grad = [];
average_grad= [];
trend = [];
pVec = [];
AverageStep = 500;
for i = 1:NTemp
    logm_smooth = smooth(log(meanMagnetization(t_range,i)), 100);
    logm_unsmooth = log(meanMagnetization(t_range,i));
    
    p = polyfit(logt, logm_unsmooth, 2);
    pVec(i) = p(1);
    
    
end

pVec

%%---------------------------------------------------------------------------------
%% qubic curvature errorbar
% fig = figure;
% fig.Position  = [0 0 300 300];
t_max = max_count;
t_range = 100:1000;
Nboot = 500;
NTemp = length(Temperature);
pVec = zeros(NTemp, Nboot);
for iboot = 1:Nboot
    BootBins = randi(Nbin, Nbin, 1);
    BootBins = BootBins';
    logt = log(t_range);
    logt = logt';
    for i = 1:NTemp
        mboot = zeros(size(squeeze(mtime1Bins(t_range, i, 1))));
        for Bi = BootBins
            mboot = mboot + squeeze(mtime1Bins(t_range, i, Bi));
        end
        mboot = mboot / Nbin;
        logm = log(mboot);
        p = polyfit(logt, logm, 2);
        pVec(i, iboot) = p(1); 
    end
end

pVecMean = mean(pVec,2);
pVecError = std(pVec, 0, 2);
yyaxis right
errorbar(Temperature, pVecMean, pVecError, ...
    'Marker', '.','MarkerSize', 20, 'LineWidth', 0.6, 'Color', 'r'); 

% line(xlim(), [0,0] ,'LineWidth', 0.5, 'Color', 'k');
grid on;
axis padded
xlabel '$T$'; 
ylabel 'Coefficient of quadratic fit'
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
% cmap = linspecer(length(Temperature));
% ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
% saveas(gcf, [save_address '/', 'Coefficient of quadratic fit.png']);
exportgraphics(gcf,[save_address '/', 'Coefficient of quadratic fit.png'],'Resolution',1200)

%%---------------------------------------------------------------------------------
%% interpol
% Assume x1vec, x2vec, y1vec, y2vec are defined as columns vectors of the same length
T1index = 4;
T2index = 6;
y1vec = log(meanMagnetization(:, T1index));
y2vec = log(meanMagnetization(:, T2index));

TPoints = 100;
n = length(y1vec); % number of interpolation intervals
x1 = Temperature(T1index);
x2 = Temperature(T2index);

% Initialize a matrix to hold the interpolated y values for each of the 100 points in each interval
interpolated_y_values = zeros(n, TPoints);
x_values = linspace(x1, x2, TPoints);
% Loop over each pair of points
for i = 1:n
    % Generate 100 equally spaced points between x1vec(i) and x2vec(i)
    
    
    % Perform linear interpolation between the two points
    interpolated_y_values(i, :) = interp1([x1, x2], [y1vec(i), y2vec(i)], x_values);
end


t_max = max_count;
t_range = 70:t_max;
logt = log(t_range);
logt = logt';
logm = interpolated_y_values(t_range,:);
% logt_smooth = smooth(log(t_range), 100);

R2Vec = zeros(size(x_values));
slope = zeros(size(x_values));
slopeErr = zeros(size(x_values));
for i = 1:length(x_values)
    % Fit the linear model
    A = logm(:,i);
    mdl = fitlm(logt, A , 'linear');
    
    % Get the R-squared value
    r_squared = mdl.Rsquared.Ordinary;
    R2Vec(i) = r_squared;
    slope(i) = mdl.Coefficients.Estimate(2);
    slopeErr(i) = mdl.Coefficients.SE(2);
end

[~, idx] = max(R2Vec);
% errorbar(Temperature, slope, slopeErr, '.','MarkerSize', 15 ); xlabel '$T$'; ylabel '$slope$'
plot(x_values, R2Vec, 'Marker', '.','MarkerSize', 15 ); xlabel '$T$'; ylabel '$R^2$'
% plot(SystemLength, R2Vec, 'Marker', '.','MarkerSize', 15 ); xlabel '$L$'; ylabel '$R^2$'

Tc = x_values(idx)
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
% cmap = linspecer(length(Temperature));
% ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'R2 interpolated.png']);


%%---------------------------------------------------------------------------------
%% interpol errorbar
% Assume x1vec, x2vec, y1vec, y2vec are defined as columns vectors of the same length
T1index = 2;
T2index = 3;

TcBin = [];
for Bi = 1:5
y1vec = log(squeeze(mtime1Bins(:, T1index, Bi)));
y2vec = log(squeeze(mtime1Bins(:, T2index, Bi)));

TPoints = 100;
n = length(y1vec); % number of interpolation intervals
x1 = Temperature(T1index);
x2 = Temperature(T2index);

% Initialize a matrix to hold the interpolated y values for each of the 100 points in each interval
interpolated_y_values = zeros(n, TPoints);
x_values = linspace(x1, x2, TPoints);
% Loop over each pair of points
for i = 1:n
    % Generate 100 equally spaced points between x1vec(i) and x2vec(i)
    
    
    % Perform linear interpolation between the two points
    interpolated_y_values(i, :) = interp1([x1, x2], [y1vec(i), y2vec(i)], x_values);
end


t_max = 1000;
t_range = 70:t_max;
logt = log(t_range);
logt = logt';
logm = (interpolated_y_values(t_range,:));
% logt_smooth = smooth(log(t_range), 100);

R2Vec = zeros(size(x_values));
slope = zeros(size(x_values));
slopeErr = zeros(size(x_values));
for i = 1:length(x_values)
    % Fit the linear model
    A = logm(:,i);
    mdl = fitlm(logt, A , 'linear');
    
    % Get the R-squared value
    r_squared = mdl.Rsquared.Ordinary;
    R2Vec(i) = r_squared;
    slope(i) = mdl.Coefficients.Estimate(2);
    slopeErr(i) = mdl.Coefficients.SE(2);
end

[~, idx] = max(R2Vec);
% errorbar(Temperature, slope, slopeErr, '.','MarkerSize', 15 ); xlabel '$T$'; ylabel '$slope$'
% plot(x_values, R2Vec, 'Marker', '.','MarkerSize', 15 ); xlabel '$T$'; ylabel '$R^2$'
% plot(SystemLength, R2Vec, 'Marker', '.','MarkerSize', 15 ); xlabel '$L$'; ylabel '$R^2$'

TcBin(Bi) = x_values(idx);
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
% cmap = linspecer(length(Temperature));
% ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
end

TcMean = mean(TcBin)
TcError = std(TcBin) / sqrt(length(TcBin))

%%---------------------------------------------------------------------------------
%% interpol spline
% Assume x1vec, x2vec, y1vec, y2vec are defined as columns vectors of the same length
T1index = 1;
T2index = 2;
T3index = 3;
T4index = 4;
y1vec = meanMagnetization(:, T1index);
y2vec = meanMagnetization(:, T2index);
y3vec = meanMagnetization(:, T3index);
y4vec = meanMagnetization(:, T4index);

TPoints = 100;
n = length(y1vec); % number of interpolation intervals
x1 = Temperature(T1index);
x2 = Temperature(T2index);
x3 = Temperature(T3index);
x4 = Temperature(T4index);

% Initialize a matrix to hold the interpolated y values for each of the 100 points in each interval
interpolated_y_values = zeros(n, TPoints);
x_values = linspace(x1, x4, TPoints);
% Loop over each pair of points
for i = 1:n
    % Generate 100 equally spaced points between x1vec(i) and x2vec(i)
    xis = [x1, x2, x3, x4];
    yis = [y1vec(i), y2vec(i), y3vec(i), y4vec(i)];
    % Perform linear interpolation between the two points
    interpolated_y_values(i, :) = interp1(xis, yis, x_values, 'spline');
end


t_max = 1000;
t_range = 40:t_max;
logt = log(t_range);
logt = logt';
logm = log(interpolated_y_values(t_range,:));
% logt_smooth = smooth(log(t_range), 100);

R2Vec = zeros(size(x_values));
slope = zeros(size(x_values));
slopeErr = zeros(size(x_values));
for i = 1:length(x_values)
    % Fit the linear model
    A = logm(:,i);
    mdl = fitlm(logt, A , 'linear');
    
    % Get the R-squared value
    r_squared = mdl.Rsquared.Ordinary;
    R2Vec(i) = r_squared;
    slope(i) = mdl.Coefficients.Estimate(2);
    slopeErr(i) = mdl.Coefficients.SE(2);
end

[~, idx] = max(R2Vec);
% errorbar(Temperature, slope, slopeErr, '.','MarkerSize', 15 ); xlabel '$T$'; ylabel '$slope$'
plot(x_values, R2Vec, 'Marker', '.','MarkerSize', 15 ); xlabel '$T$'; ylabel '$R^2$'
% plot(SystemLength, R2Vec, 'Marker', '.','MarkerSize', 15 ); xlabel '$L$'; ylabel '$R^2$'

Tc = x_values(idx)
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
% cmap = linspecer(length(Temperature));
% ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

%%---------------------------------------------------------------------------------
%% interpol spline errorbar
% Assume x1vec, x2vec, y1vec, y2vec are defined as columns vectors of the same length
T1index = 1;
T2index = 2;
T3index = 3;
T4index = 4;

TcBin = [];
for Bi = 1:5
y1vec = squeeze(mtime1Bins(:, T1index, Bi));
y2vec = squeeze(mtime1Bins(:, T2index, Bi));
y3vec = squeeze(mtime1Bins(:, T3index, Bi));
y4vec = squeeze(mtime1Bins(:, T4index, Bi));

TPoints = 100;
n = length(y1vec); % number of interpolation intervals
x1 = Temperature(T1index);
x2 = Temperature(T2index);
x3 = Temperature(T3index);
x4 = Temperature(T4index);

% Initialize a matrix to hold the interpolated y values for each of the 100 points in each interval
interpolated_y_values = zeros(n, TPoints);
x_values = linspace(x1, x4, TPoints);
% Loop over each pair of points
for i = 1:n
    % Generate 100 equally spaced points between x1vec(i) and x2vec(i)
    
    xis = [x1, x2, x3, x4];
    yis = [y1vec(i), y2vec(i), y3vec(i), y4vec(i)];
    % Perform linear interpolation between the two points
    interpolated_y_values(i, :) = interp1(xis, yis, x_values, 'spline');
end


t_max = max_count;
t_range = 100:t_max;
logt = log(t_range);
logt = logt';
logm = log(interpolated_y_values(t_range,:));
% logt_smooth = smooth(log(t_range), 100);

R2Vec = zeros(size(x_values));
slope = zeros(size(x_values));
slopeErr = zeros(size(x_values));
for i = 1:length(x_values)
    % Fit the linear model
    A = logm(:,i);
    mdl = fitlm(logt, A , 'linear');
    
    % Get the R-squared value
    r_squared = mdl.Rsquared.Ordinary;
    R2Vec(i) = r_squared;
    slope(i) = mdl.Coefficients.Estimate(2);
    slopeErr(i) = mdl.Coefficients.SE(2);
end

[~, idx] = max(R2Vec);
% errorbar(Temperature, slope, slopeErr, '.','MarkerSize', 15 ); xlabel '$T$'; ylabel '$slope$'
% plot(x_values, R2Vec, 'Marker', '.','MarkerSize', 15 ); xlabel '$T$'; ylabel '$R^2$'
% plot(SystemLength, R2Vec, 'Marker', '.','MarkerSize', 15 ); xlabel '$L$'; ylabel '$R^2$'

TcBin(Bi) = x_values(idx);
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
% cmap = linspecer(length(Temperature));
% ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
end

TcMean = mean(TcBin)
TcError = std(TcBin) / sqrt(length(TcBin))


%%---------------------------------------------------------------------------------
%% segments slope
t_max = max_count;
t_range = 70:t_max;

x = t_range;
y = meanMagnetization(t_range,1);

% Generate data
% x = linspace(1, 1000, 1000); % 1000 points from 1 to 1000
% y = x.^2 + 10 * (1+randn(size(x))); % Example y data with some noise


% Log-transform the data
log_x = log(x)';
log_y = log(y)';

% Number of segments
n_segments = 2;

% Generate log-spaced indices
start_segmentation_index = 1;
segment_boundaries = linspace(log_x(start_segmentation_index), log_x(end), n_segments+1);

% Initialize array for slopes
slopes = zeros(1, n_segments);
midpoints = [];
% Calculate slopes for each segment
for i = 1:n_segments
    % Find indices for the current segment
    segment_start = find(log_x >= segment_boundaries(i), 1, 'first');
    segment_end = find(log_x <= segment_boundaries(i+1), 1, 'last');
    
    x_segment = log_x(segment_start:segment_end);
    y_segment = log_y(segment_start:segment_end);
    
    % Perform linear fit (1st degree polynomial) on log-log data
    p = polyfit(x_segment, y_segment, 1);
    slopes(i) = p(1); % Slope of the segment

    % midpoints(i) = (segment_boundaries(i) + segment_boundaries(i+1)) / 2;
    midpoints(i) = (exp(segment_boundaries(i)) + exp(segment_boundaries(i+1))) / 2;

end

% Calculate the variance of the slopes
slope_variance = var(slopes);

% Display the results
disp([sprintf('Variance of the slopes: \n'), num2str(slope_variance)]);
% disp([sprintf('Slopes: \n'), num2str(slopes)]);

% Plot the slopes as a function of 1/midpoints
inverse_midpoints = 1 ./ midpoints;
figure;
plot(inverse_midpoints, slopes, 'o-');
xlabel('1 / Midpoint of Segment');
ylabel('Slope');
title('Slope vs. 1/Midpoint of Segment');
grid on;

% % Plot the original data and the log-log segments
% figure;
% loglog(x(start_segmentation_index:end), y(start_segmentation_index:end), 'b-'); % Original data on log-log scale
% hold on;
% for i = 1:n_segments
%     segment_start = find(log_x >= segment_boundaries(i), 1, 'first');
%     segment_end = find(log_x <= segment_boundaries(i+1), 1, 'last');
% 
%     x_segment = log_x(segment_start:segment_end);
%     y_segment = log_y(segment_start:segment_end);
% 
%     % Plot segments on log-log scale
%     plot(exp(x_segment), exp(y_segment), 'r-', 'LineWidth', 1.5);
% end
% hold off;
% xlabel('x');
% ylabel('y');
% title('Log-Log Plot with Log-Spaced Segments');
% legend('Original Data', 'Segmented Fits');

%%---------------------------------------------------------------------------------
%% segments slope multiple
t_max = max_count;
t_range = 70:1000;

NSystem = length(Temperature);
for k = 1:NSystem

x = t_range;
y = meanMagnetization(t_range,k);

% Log-transform the data
log_x = log(x)';
log_y = log(y)';

% Number of segments
n_segments = 2;

% Generate log-spaced indices
start_segmentation_index = 1;
segment_boundaries = linspace(log_x(start_segmentation_index), log_x(end), n_segments+1);

% Initialize array for slopes
slopes = zeros(1, n_segments);
midpoints = [];
% Calculate slopes for each segment
for i = 1:n_segments
    % Find indices for the current segment
    segment_start = find(log_x >= segment_boundaries(i), 1, 'first');
    segment_end = find(log_x <= segment_boundaries(i+1), 1, 'last');
    
    x_segment = log_x(segment_start:segment_end);
    y_segment = log_y(segment_start:segment_end);
    
    % Perform linear fit (1st degree polynomial) on log-log data
    p = polyfit(x_segment, y_segment, 1);
    slopes(i) = p(1); % Slope of the segment

    % midpoints(i) = (segment_boundaries(i) + segment_boundaries(i+1)) / 2;
    midpoints(i) = (exp(segment_boundaries(i)) + exp(segment_boundaries(i+1))) / 2;

end

% Calculate the variance of the slopes
slope_variance = var(slopes);

% Plot the slopes as a function of 1/midpoints
inverse_midpoints = 1 ./ midpoints;
plot(inverse_midpoints, slopes, 'o-');
hold on;
end

xlabel('1 / Midpoint of Segment');
ylabel('Segment slope');
% a = xlim;
% xlim([0 a(2)]);
legendStrings = "$L$ = " + string(SystemLength) + ", $T$ = " + string(Temperature);
legend(legendStrings, 'Location','best', ...
    'Interpreter','latex');
grid on;
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));

ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'Segment slope with ' num2str(n_segments) ...
            ' segments.png']);

%%---------------------------------------------------------------------------------
%% interpol segment slpoe
% Assume x1vec, x2vec, y1vec, y2vec are defined as columns vectors of the same length
T1index = 2;
T2index = 3;
y1vec = log(meanMagnetization(:, T1index));
y2vec = log(meanMagnetization(:, T2index));

TPoints = 100;
n = length(y1vec); % number of interpolation intervals
x1 = Temperature(T1index);
x2 = Temperature(T2index);

% Initialize a matrix to hold the interpolated y values for each of the 100 points in each interval
interpolated_y_values = zeros(n, TPoints);
x_values = linspace(x1, x2, TPoints);

% Loop over each pair of points
for i = 1:n
    % Generate 100 equally spaced points between y1vec(i) and y2vec(i)
    % Perform linear interpolation between the two points
    interpolated_y_values(i, :) = interp1([x1, x2], [y1vec(i), y2vec(i)], x_values);
    
    
end

t_max = max_count;
t_range = 70:t_max;
logt = log(t_range);
logt = logt';
logm = interpolated_y_values(t_range,:);
% logt_smooth = smooth(log(t_range), 100);

R2Vec = zeros(size(x_values));
slope = zeros(size(x_values));
slopeErr = zeros(size(x_values));
for i = 1:length(x_values)
    % Log-transform the data
    log_x = logt;
    log_y = logm(:,i);

    % Number of segments
    n_segments = 2;

    % Generate log-spaced indices
    start_segmentation_index = 1;
    segment_boundaries = linspace(log_x(start_segmentation_index), log_x(end), n_segments+1);

    % Initialize array for slopes
    segmentSlopes = zeros(1, n_segments);
    midpoints = [];
    % Calculate slopes for each segment
    for j = 1:n_segments
        % Find indices for the current segment
        segment_start = find(log_x >= segment_boundaries(j), 1, 'first');
        segment_end = find(log_x <= segment_boundaries(j+1), 1, 'last');

        x_segment = log_x(segment_start:segment_end);
        y_segment = log_y(segment_start:segment_end);

        % Perform linear fit (1st degree polynomial) on log-log data
        p = polyfit(x_segment, y_segment, 1);
        segmentSlopes(j) = p(1); % Slope of the segment

        % midpoints(i) = (segment_boundaries(i) + segment_boundaries(i+1)) / 2;
        midpoints(j) = (exp(segment_boundaries(j)) + exp(segment_boundaries(j+1))) / 2;

    end

    % Calculate the variance of the slopes
    segment_slope_variance = var(segmentSlopes);

    % Plot the slopes as a function of 1/midpoints
    inverse_midpoints = 1 ./ midpoints;


    % Fit the linear model
    mdl = fitlm(inverse_midpoints, segmentSlopes , 'linear');

    % Get the R-squared value
    r_squared = mdl.Rsquared.Ordinary;
    R2Vec(i) = r_squared;
    slope(i) = mdl.Coefficients.Estimate(2);
    slopeErr(i) = mdl.Coefficients.SE(2);
end

[~, idx] = min(abs(slope));
% errorbar(x_values, (slope), slopeErr, '.','MarkerSize', 15 ); xlabel '$T$'; ylabel 'trend of segment slopes'
errorbar(x_values, abs(slope), slopeErr, '.','MarkerSize', 15 ); xlabel '$T$'; ylabel 'trend of segment slopes'
% plot(x_values, R2Vec, 'Marker', '.','MarkerSize', 15 ); xlabel '$T$'; ylabel 'trend of segment slopes'
% plot(SystemLength, R2Vec, 'Marker', '.','MarkerSize', 15 ); xlabel '$L$'; ylabel '$R^2$'
grid on;
Tc = x_values(idx)
set(gca,'fontsize',12)
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
% cmap = linspecer(length(Temperature));
% ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'trend of segment slopes interpolated.png']);



%%---------------------------------------------------------------------------------
%% overlapping segments slope
t_max = max_count;
t_range = 80:t_max;

x = t_range;
y = meanMagnetization(t_range, 2);

% % Generate data
% x = linspace(1, 1000, 1000); % 1000 points from 1 to 1000
% y = x.^2 + 10 * (1+randn(size(x))); % Example y data with some noise

% Log-transform the data
log_x = log(x)';
log_y = log(y)';

% Number of segments
n_segments = 20;

% Calculate the total length in log scale
total_length = log_x(end) - log_x(1);

% Calculate the length of each segment
segment_length = total_length / (3); % Each segment covers half of the total length

% Calculate the step size between segment starts
step_size = (total_length - segment_length) / (n_segments - 1);

% Initialize arrays for slopes and midpoints
slopes = zeros(1, n_segments);
midpoints = zeros(1, n_segments);

% Calculate slopes for each overlapping segment
for i = 1:n_segments
    % Calculate start and end of the current segment
    % segment_start = find(log_x >= segment_boundaries(i), 1, 'first');
    % segment_end = find(log_x <= segment_boundaries(i+1), 1, 'last');
    % % 
    % x_segment = log_x(segment_start:segment_end);
    % y_segment = log_y(segment_start:segment_end);

    % 
    segment_start = log_x(1) + (i-1) * step_size;
    segment_end = segment_start + segment_length;
    
    assert(segment_start >= log_x(1));
    assert(segment_end <=log_x(end));
    
    % Find indices for the current segment
    start_index = find(log_x >= segment_start, 1, 'first');
    end_index = find(log_x <= segment_end, 1, 'last');

    x_segment = log_x(start_index:end_index);
    y_segment = log_y(start_index:end_index);
    
    % Perform linear fit (1st degree polynomial) on log-log data
    p = polyfit(x_segment, y_segment, 1);
    slopes(i) = p(1); % Slope of the segment

    % Calculate midpoint in original scale
    midpoints(i) = exp((segment_start + segment_end) / 2);
end

% Calculate the variance of the slopes
slope_variance = var(slopes);

% Display the results
disp([sprintf('Variance of the slopes: \n'), num2str(slope_variance)]);

% Plot the slopes as a function of 1/midpoints
inverse_midpoints = 1 ./ midpoints;
figure;
plot(inverse_midpoints, slopes, 'o-');
xlabel('1 / Midpoint of Segment');
ylabel('Slope');
title('Slope vs. 1/Midpoint of Segment');
grid on;

% Optional: Visualize the segments
% figure;
% loglog(x, y, 'b.');
% hold on;
% for i = 1:n_segments
%     segment_start = log_x(1) + (i-1) * step_size;
%     segment_end = segment_start + segment_length;
%     plot([segment_start, segment_end], [], 'r-', 'LineWidth', 2);
% end
% xlabel('x');
% ylabel('y');
% title('Data with Overlapping Segments');
% hold off;

%%---------------------------------------------------------------------------------
%% overlapping segments slope multiple
t_max = max_count;
t_range = 100:1000;

NSystem = length(Temperature);
for k = 1:NSystem

    x = t_range;
    y = meanMagnetization(t_range,k);

    % Log-transform the data
    log_x = log(x)';
    log_y = log(y)';

    % Number of segments
    n_segments = 40;

    % Calculate the total length in log scale
    total_length = log_x(end) - log_x(1);

    % Calculate the length of each segment
    segment_length = total_length / (2); % Each segment covers half of the total length

    % Calculate the step size between segment starts
    step_size = (total_length - segment_length) / (n_segments - 1);

    % Initialize arrays for slopes and midpoints
    slopes = zeros(1, n_segments);
    slopesErr = zeros(1, n_segments);
    midpoints = zeros(1, n_segments);

    % Calculate slopes for each overlapping segment
    for i = 1:n_segments
        % Calculate start and end of the current segment
        % segment_start = find(log_x >= segment_boundaries(i), 1, 'first');
        % segment_end = find(log_x <= segment_boundaries(i+1), 1, 'last');
        % %
        % x_segment = log_x(segment_start:segment_end);
        % y_segment = log_y(segment_start:segment_end);

        %
        segment_start = log_x(1) + (i-1) * step_size;
        segment_end = segment_start + segment_length;

        assert(segment_start >= log_x(1));
        assert(segment_end <=log_x(end));

        % Find indices for the current segment
        start_index = find(log_x >= segment_start, 1, 'first');
        end_index = find(log_x <= segment_end, 1, 'last');

        x_segment = log_x(start_index:end_index);
        y_segment = log_y(start_index:end_index);

        % Perform linear fit (1st degree polynomial) on log-log data
        % p = polyfit(x_segment, y_segment, 1);
        % slopes(i) = p(1); % Slope of the segment

        mdl = fitlm(x_segment, y_segment , 'linear');
        % lambda_m = mdl.Coefficients.Estimate(2);

        slopes(i) = mdl.Coefficients.Estimate(2); % Slope of the segment
        slopesErr(i) = mdl.Coefficients.SE(2); % Slope of the segment

        % Calculate midpoint in original scale
        midpoints(i) = exp((segment_start + segment_end) / 2);
    end
    assert(end_index==length(t_range));
    % Calculate the variance of the slopes
    slope_variance = var(slopes);

    % Plot the slopes as a function of 1/midpoints
    inverse_midpoints = 1 ./ midpoints;
    % plot(inverse_midpoints, slopes, 'o-');
    errorbar(inverse_midpoints, slopes, slopesErr, 'o-');

    hold on;
end

xlabel('1 / Midpoint of Segment');
ylabel('Segment slope');
% a = xlim;
% xlim([0 a(2)]);
legendStrings = "$L$ = " + string(SystemLength) + ", $T$ = " + string(Temperature);
% legendStrings = "$T$ = " + string(Temperature);
legend(legendStrings, 'Location','best', ...
    'Interpreter','latex');
grid on;
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));
% set(gca, 'YScale', 'log')   %using plot


ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
% saveas(gcf, [save_address '/', 'running segment slope.png']);

saveas(gcf, [save_address '/', 'Segment slope with ' num2str(n_segments) ...
            ' overlapping segments and length ' num2str(segment_length/ total_length) ' of total logarithmic length .png']);

%%---------------------------------------------------------------------------------
%% overlapping segments slope multiple error bar
t_max = max_count;
t_range = 100:1000;
NSystem = length(Temperature);
Nboot = 100;
for k = 1:NSystem

    x = t_range;
    % y = meanMagnetization(t_range,k);

    % Log-transform the data
    log_x = log(x)';
    % log_y = log(y)';

    % Number of segments
    n_segments = 7;

    % Calculate the total length in log scale
    total_length = log_x(end) - log_x(1);

    % Calculate the length of each segment
    segment_length = total_length / (2); % Each segment covers half of the total length

    % Calculate the step size between segment starts
    step_size = (total_length - segment_length) / (n_segments - 1);

    % Initialize arrays for slopes and midpoints
    slopes = zeros(1, n_segments);
    slopesErr = zeros(1, n_segments);
    midpoints = zeros(1, n_segments);

    % Calculate slopes for each overlapping segment
    for i = 1:n_segments
        % Calculate start and end of the current segment
        % segment_start = find(log_x >= segment_boundaries(i), 1, 'first');
        % segment_end = find(log_x <= segment_boundaries(i+1), 1, 'last');
        % %
        % x_segment = log_x(segment_start:segment_end);
        % y_segment = log_y(segment_start:segment_end);

        %
        segment_start = log_x(1) + (i-1) * step_size;
        segment_end = segment_start + segment_length;

        assert(segment_start >= log_x(1));
        assert(segment_end <=log_x(end));

        % Find indices for the current segment
        start_index = find(log_x >= segment_start, 1, 'first');
        end_index = find(log_x <= segment_end, 1, 'last');

        x_segment = log_x(start_index:end_index);
        % y_segment = log_y(start_index:end_index);

        lambda_mBin = [];
        for iboot = 1:Nboot
            BootBins = randi(Nbin, Nbin, 1);
            BootBins = BootBins';
            mboot = zeros(size(squeeze(mtime1Bins(t_range, k, 1))));
            for Bi = BootBins
                mboot = mboot + squeeze(mtime1Bins(t_range, k, Bi));
            end
            mboot = mboot / Nbin;
            % logm = log(mboot);
            mdl = fitlm(log_x(start_index:end_index), log(mboot(start_index:end_index)) , 'linear');
            lambda_m = mdl.Coefficients.Estimate(2);

            lambda_mBin(iboot) = lambda_m;
        end
        lambda_m_mean = mean(lambda_mBin);
        lambda_m_error = std(lambda_mBin);

        %
        % disp(['z = ' fmtMeanUnc(z_mean, z_error)]);
        % disp(['nu = ' fmtMeanUnc(nu_mean, nu_error)]);
        % disp(['beta = ' fmtMeanUnc(beta_mean, beta_error)]);
        % disp(['alpha = ' computeAndFormatMean(alphaBin)]);
        % % disp(['eta = ' computeAndFormatMean(etaBin), newline]);

        slopes(i) = lambda_m_mean; % Slope of the segment
        slopesErr(i) = lambda_m_error;

        % Calculate midpoint in original scale
        midpoints(i) = exp((segment_start + segment_end) / 2);
    end
    assert(end_index==length(t_range));
    % Calculate the variance of the slopes
    slope_variance = var(slopes);

    % Plot the slopes as a function of 1/midpoints
    % inverse_midpoints = 1 ./ midpoints;
    inverse_midpoints = midpoints;
    errorbar(inverse_midpoints, slopes, slopesErr, 'o-');
    hold on;
end
xlabel('Midpoint of Segment');
% xlabel('1 / Midpoint of Segment');
ylabel('Segment slope');
% a = xlim;
ylim([-0.090 -0.074]);
legendStrings = "$L$ = " + string(SystemLength) + ", $T$ = " + string(Temperature);
% legendStrings = "$T$ = " + string(Temperature);
legend(legendStrings, 'Location','northwest', ...
    'Interpreter','latex');
grid on;
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));
% set(gca, 'YScale', 'log')   %using plot


ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
% saveas(gcf, [save_address '/', 'running segment slope.png']);

saveas(gcf, [save_address '/', 'Segment slope with ' num2str(n_segments) ...
            ' overlapping segments and length ' num2str(segment_length/ total_length) ' of total logarithmic length .png']);


%%---------------------------------------------------------------------------------
%% overlapping segments critical exponents
t_max = max_count;
t_range = 70:t_max;
i_calc = 3;

x = t_range;
y = meanMagnetization(t_range, 1);

% Log-transform the data
log_x = log(x)';
log_y = log(y)';

% Number of segments
n_segments = 20;

% Calculate the total length in log scale
total_length = log_x(end) - log_x(1);

% Calculate the length of each segment
segment_length = total_length / (1.5); % Each segment covers half of the total length

% Calculate the step size between segment starts
step_size = (total_length - segment_length) / (n_segments - 1);

% Initialize arrays for slopes and midpoints
slopes = zeros(1, n_segments);
midpoints = zeros(1, n_segments);
slopesErr = zeros(1, n_segments);
% Calculate slopes for each overlapping segment
for i = 1:n_segments
    % Calculate start and end of the current segment
    segment_start = log_x(1) + (i-1) * step_size;
    segment_end = segment_start + segment_length;

    assert(segment_start >= log_x(1));
    assert(segment_end <=log_x(end));

    % Find indices for the current segment
    start_index = find(log_x >= segment_start, 1, 'first');
    end_index = find(log_x <= segment_end, 1, 'last');

    x_segment = log_x(start_index:end_index);
    y_segment = log_y(start_index:end_index);

    logt = log(t_range);
    logt = logt';

    zBin = [];
    nuBin = [];
    betaBin = [];
    alphaBin = [];

    lambda_mBin = [];
    lambda_mmBin = [];
    lambda_meBin = [];
    lambda_eeBin = [];
    for Bi = 1:5
        logt = log(t_range);
        logt = logt';
        m1 = squeeze(mtime1Bins(t_range, i_calc, Bi));
        m2 = squeeze(mtime2Bins(t_range, i_calc, Bi));
        m3 = squeeze(mtime3Bins(t_range, i_calc, Bi));
        m4 = squeeze(mtime4Bins(t_range, i_calc, Bi));
        E1 = squeeze(Etime1Bins(t_range, i_calc, Bi));
        E2 = squeeze(Etime2Bins(t_range, i_calc, Bi));
        E3 = squeeze(Etime3Bins(t_range, i_calc, Bi));
        E4 = squeeze(Etime4Bins(t_range, i_calc, Bi));
        m1E1 = squeeze(m1E1timeBins(t_range, i_calc, Bi));

        fmmBin = m2 ./ (m1 .^2) - 1;
        feeBin = E2 ./ (E1 .^2) - 1;
        fmeBin = m1E1 ./ (m1 .* E1) - 1;

        % p = polyfit(logt, log(m1), 1);
        % lambda_m = p(1);
        %
        % p = polyfit(logt, log(fmmBin), 1);
        % lambda_mm = p(1);
        %
        % p = polyfit(logt, log(fmeBin), 1);
        % lambda_me = p(1);
        %
        % p = polyfit(logt, log(feeBin), 1);
        % lambda_ee = p(1);

        mdl = fitlm(logt(start_index:end_index), log(m1(start_index:end_index)) , 'linear');
        lambda_m = mdl.Coefficients.Estimate(2);

        mdl = fitlm(logt(start_index:end_index), log(fmmBin(start_index:end_index)) , 'linear');
        lambda_mm = mdl.Coefficients.Estimate(2);

        mdl = fitlm(logt(start_index:end_index), log(fmeBin(start_index:end_index)) , 'linear');
        lambda_me = mdl.Coefficients.Estimate(2);

        mdl = fitlm(logt(start_index:end_index), log(feeBin(start_index:end_index)) , 'linear');
        lambda_ee = mdl.Coefficients.Estimate(2);

        zBin(Bi) = 2 / lambda_mm;
        nuBin(Bi) = lambda_mm / (2 * lambda_me);
        betaBin(Bi) = lambda_m / lambda_me;
        alphaBin(Bi) = lambda_ee / lambda_me;


        lambda_mBin(Bi) = lambda_m;
        lambda_mmBin(Bi) = lambda_mm;
        lambda_meBin(Bi) = lambda_me;
        lambda_eeBin(Bi) = lambda_ee;
    end
    z_mean = mean(zBin);
    nu_mean = mean(nuBin);
    beta_mean = mean(betaBin);
    alphamean = mean(alphaBin);  
    z_error = std(zBin, 0) / sqrt(length(zBin));
    nu_error = std(nuBin, 0) / sqrt(length(nuBin));
    beta_error = std(betaBin, 0) / sqrt(length(betaBin));
    % 
    % disp(['z = ' fmtMeanUnc(z_mean, z_error)]);
    % disp(['nu = ' fmtMeanUnc(nu_mean, nu_error)]);
    % disp(['beta = ' fmtMeanUnc(beta_mean, beta_error)]);
    % disp(['alpha = ' computeAndFormatMean(alphaBin)]);
    % % disp(['eta = ' computeAndFormatMean(etaBin), newline]);

    % Perform linear fit (1st degree polynomial) on log-log data
    p = polyfit(x_segment, y_segment, 1);
    slopes(i) = nu_mean; % Slope of the segment
    slopesErr(i) = nu_error;
    % Calculate midpoint in original scale
    midpoints(i) = exp((segment_start + segment_end) / 2);
end

% Calculate the variance of the slopes
slope_variance = var(slopes);

% Display the results
disp([sprintf('Variance of the slopes: \n'), num2str(slope_variance)]);

% Plot the slopes as a function of 1/midpoints
inverse_midpoints = 1 ./ midpoints;
figure;
% plot(inverse_midpoints, slopes, 'o-');
errorbar(inverse_midpoints, slopes, slopesErr, 'o-');

xlabel('1 / Midpoint of Segment');
ylabel('Slope');
title('Slope vs. 1/Midpoint of Segment');
grid on;

p = polyfit(inverse_midpoints, slopes, 1);
stdSlopes = std(slopes)
intercept_at_x0 = polyval(p, 0)

saveas(gcf, [save_address '/', 'beta.png']);

% Optional: Visualize the segments
% figure;
% loglog(x, y, 'b.');
% hold on;
% for i = 1:n_segments
%     segment_start = log_x(1) + (i-1) * step_size;
%     segment_end = segment_start + segment_length;
%     plot([segment_start, segment_end], [], 'r-', 'LineWidth', 2);
% end
% xlabel('x');
% ylabel('y');
% title('Data with Overlapping Segments');
% hold off;

%%---------------------------------------------------------------------------------
%% data collapse
clf('reset')
t_max = length(meanMagnetization);
t_range = 1:t_max;

for i = 1:length(Temperature)
    % errorbar(t_range, meanMagnetization(t_range,i), errorMagnetization(t_range,i), ...
    %     'Marker', '.', 'MarkerSize', 5);
    XX = t_range.^(0.549).*(Temperature(i)/0.8503 - 1);
    YY = meanMagnetization(t_range,i).*t_range.^(-0.0837);
    % plot(t_range, errorMagnetization(t_range,i), 'Marker', '.', 'MarkerSize', 7);
    plot(XX, YY, 'Marker', '.', 'MarkerSize', 7);
    % plot(log(t_range), log(meanMagnetization(t_range,i)), 'Marker', '.', 'MarkerSize', 7);

    hold on
end


xlabel '$t$'
ylabel '$|\mathbf{m}|$'
% ylabel '$\mathrm{error}(|\mathbf{m}|)$'

legendStrings = "$L$ = " + string(SystemLength) + ", $T$ = " + string(Temperature);
% legendStrings = "$T$ = " + string(Temperature);
legend(legendStrings, 'Location','best', ...
    'Interpreter','latex')
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
ax = gca;
ax.TickLength(1) = 0.02;
cmap = linspecer(length(Temperature));
ax.ColorOrder = cmap;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
% saveas(gcf, [save_address '/', 'm.png']);
% saveas(gcf, [save_address '/', 'm_scaled.png']);
% saveas(gcf, [save_address '/', 'merror.png']);
% saveas(gcf, [save_address '/', 'merror_scaled.png']);


% Tc = 0.844;
% reducedTemp = data(:,1,:) / Tc -1;
% L = sqrt(listN);
% for i=1:N
%     Xtilde(:,i) = data(:,7,i);  %Binder
%     x(:,i) = L(i)^(1) * reducedTemp(:,:,i);
% end
% 
% clf('reset')
% for i = 1:N
%     plot(x(:,i), Xtilde(:,i) , '.','MarkerSize', 12);
% %     plot(x(:,i), data(:,7,i), '.','MarkerSize', 10);
%     hold on
% end
% % set(gca,'TickLabelInterpreter','latex');
% % xlim([-5 5])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% 
% % ylabel '$m L^{\beta / \nu}$'
% % ylabel '$\chi L^{-\gamma / \nu}$'
% ylabel '$U_4$'
% xlabel('$t L^{1/ \nu}$')
% legendStrings = "L = " + string(listL);
% legend(legendStrings, 'Location','southeast','Interpreter','latex');
% set(gca,'fontsize',fontSize) 
% set(gca,'XMinorTick','on','YMinorTick','on');
% ax = gca;
% cmap = turbo(N+1);
% ax.ColorOrder = cmap;
% ax.TickLength(1) = 0.02;
% set(gca,'TickLabelInterpreter','latex');
% set(0,'defaulttextinterpreter','latex');
% exportgraphics(gcf,[save_address '/', 'Fss binder.png'],'Resolution',dpi)

%% functions
function BootAveraged = BootAverage(y, BootBins)
    % [maxY, maxIndex] = max(yy
    yBoot = zeros(size(squeeze(y(:, :, 1))));
    for l = BootBins
        yBoot = yBoot + squeeze(y(:, :, l));
    end
    BootAveraged = yBoot / length(BootBins);
end
