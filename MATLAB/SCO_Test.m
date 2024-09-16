%% load data
clc
clear all
path1 = ['E:\visualize\NewAnalysis\SCO\Validity\Random\'];
path2 = ['E:\visualize\NewAnalysis\SCO\Validity\Overrelaxed\'];
methods = ["Metropolis", "at spin update",  ...
            "1 sweep", "10 sweeps",];

% i = 1;
% for method = methods
%     MethodName = sprintf('%s',method);
%     data(:,:,i) =  readmatrix([path, MethodName, '\MeanX.csv']);
%     error(:,:,i) = readmatrix([path, MethodName, '\errorX.csv']);
%     i = i + 1;
% end


date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
save_address = ['E:/visualize/output/', date_time_string , '_SCO'];
[status, msg, msgID] = mkdir(save_address)
dpi = 600;
%%---------------------------------------------------------------------------------
%% plot m

clf('reset')
cmap = linspecer(length(methods));
markerShapes = {'v', 's', '^', 'o'};
legendStrings=[];
mi=0;
path = path1;
i = 1;
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    data(:,:,i) =  readmatrix([path, MethodName, '\MeanX.csv']);
    error(:,:,i) = readmatrix([path, MethodName, '\errorX.csv']);
    i = i + 1;
    
    DataIndex = 4;

    DataToPlot = data(:,DataIndex, mi);
    errorToPlot = error(:,DataIndex, mi);
    Temprature = data(:,1, mi);
    h = errorbar(Temprature, DataToPlot, errorToPlot, ...
            'Marker', markerShapes{mi} ,'MarkerSize', 6, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    % h = plot(Temprature, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 4
         h.MarkerFaceColor = h.Color;
    % end
    hold on
    legendStrings = [legendStrings, string([ MethodName])];

end

xlim([0.2 1.25])
% ylim([0 0.26])
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =12;
set(gca,'fontsize',myFontSize) 
ylabel ('$m$', 'FontSize', myFontSize)
xlabel ('$T$', 'FontSize', myFontSize)
leg1 = legend(legendStrings, 'Location','southwest','Interpreter','latex', 'FontSize', myFontSize);
% leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', myFontSize);

set(leg1,'Box','off')

set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
% axis(ax, 'tight');
pbaspect([2 1 1])

% 
% exportgraphics(gcf,[save_address '/', 'SCO_RandomVailidity.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'SCO_RandomVailidity.eps'])


exportgraphics(gcf,[save_address '/', 'SCO_ORVailidity.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/', 'SCO_ORVailidity.eps'])

%% load data
clc
clear all
path = ['E:\visualize\NewAnalysis\SCO\Time series\Random\'];
% path = ['E:\visualize\NewAnalysis\SCO\Validity\Overrelaxed\'];
methods = ["Metropolis", "at spin update",  ...
            "1 sweep", "10 sweeps",];

% i = 1;
% for method = methods
%     MethodName = sprintf('%s',method);
%     data(:,:,i) =  readmatrix([path, MethodName, '\MeanX.csv']);
%     error(:,:,i) = readmatrix([path, MethodName, '\errorX.csv']);
%     i = i + 1;
% end


date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
save_address = ['E:/visualize/output/', date_time_string , '_SCO'];
[status, msg, msgID] = mkdir(save_address)
dpi = 600;
%%---------------------------------------------------------------------------------
%% plot C

clf('reset')
cmap = linspecer(length(methods));
markerShapes = {'v', 's', '^', 'o'};
legendStrings=[];
mi=0;
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    data_path = [path, MethodName, '\MeanX.csv'];
    for j = 1:8
        data_path = [path, MethodName, '\', num2str(j-1), '\', 'EnsembleResults.h5'];
        mtime1Bins(:, i, j) = h5read(data_path, ['/mtime1Bins/', num2str(j-1)]);
    end
    i = i + 1;
    
    DataIndex = 5; %magnetization

    DataToPlot = data(:,DataIndex, mi);
    errorToPlot = error(:,DataIndex, mi);
    Temprature = data(:,1, mi);
    h = errorbar(Temprature, DataToPlot, errorToPlot, ...
            'Marker', markerShapes{mi} ,'MarkerSize', 6, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    % h = plot(Temprature, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 4
         h.MarkerFaceColor = h.Color;
    % end
    hold on
    legendStrings = [legendStrings, string([ MethodName])];

end

% xlim([1e-3 1e1])
% ylim([0 0.26])
yl = ylim;
xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =12;
set(gca,'fontsize',myFontSize) 
ylabel ('$m$', 'FontSize', myFontSize)
xlabel ('$T$', 'FontSize', myFontSize)
leg1 = legend(legendStrings, 'Location','southwest','Interpreter','latex', 'FontSize', myFontSize);
set(leg1,'Box','off')

set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
axis(ax, 'tight');
pbaspect([2 1 1])


exportgraphics(gcf,[save_address '/', 'SCO_ORVailidity.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/', 'SCO_ORVailidity.eps'])
