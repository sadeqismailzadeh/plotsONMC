%% load data
clc
clear all
path = ['E:\visualize\NewAnalysis\DynamicThinning_Vs_Fukui\Temeprature\'];
methods = ["Fukui-Todo", "dyn. thinning"];

i = 1;
for method = methods
    MethodName = sprintf('%s',method);
    data(:,:,i) =  readmatrix([path, MethodName, '\MeanX.csv']);
    error(:,:,i) = readmatrix([path, MethodName, '\errorX.csv']);
    i = i + 1;
end


date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
save_address = ['E:/visualize/output/', date_time_string , '_DythinVsFukui'];
[status, msg, msgID] = mkdir(save_address)
dpi = 600;
%% set subfigs different OR Clock Main
fig = figure;
fig.Position  = [0 0 500 180];
% tiledlayout(4,2)
axs = tight_subplot(1,2,[.05 .13],[.25 .03],[.13 .03])

%%---------------------------------------------------------------------------------
%% plot C

% f = figure;
% f.Position = [100 100 300 200];
cmap = linspecer(length(methods));
markerShapes = {'^', 'v'};
legendStrings=[];
mi=0;
axes(axs(1))
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    data(:,:,i) =  readmatrix([path, MethodName, '\MeanX.csv']);
    error(:,:,i) = readmatrix([path, MethodName, '\errorX.csv']);
    i = i + 1;
    
    DataIndex = 31; %Pair

    DataToPlot = data(:,DataIndex, mi); %run/step
    Temprature = data(:,1, mi);
    % h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(Temprature, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 4
         h.MarkerFaceColor = h.Color;
    % end
    hold on
    legendStrings = [legendStrings, string([ MethodName])];;

end

% xlim([1e-3 1e1])
% ylim([0 0.26])
yl = ylim;
xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
xticks([1e-3, 1e-2, 1e-1, 1e0])
myFontSize =12;
set(gca,'fontsize',myFontSize) 
ylabel ('$N_\mathrm{pw,acc}$', 'FontSize', myFontSize)
xlabel ('$T$', 'FontSize', myFontSize)
leg1 = legend(legendStrings, 'Location','southwest','Interpreter','latex', 'FontSize', myFontSize);
set(leg1,'Box','off')

set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
% pbaspect([2 1 1])
set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot

% 
% exportgraphics(gcf,[save_address '/', 'DythinVsFukui_C_accept.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'DythinVsFukui_C_accept.eps'])

%%---------------------------------------------------------------------------------
%% plot RND

% f = figure;
% f.Position = [100 100 300 200];
cmap = linspecer(length(methods));
markerShapes = {'^', 'v'};
legendStrings=[];
mi=0;
axes(axs(2))
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    data(:,:,i) =  readmatrix([path, MethodName, '\MeanX.csv']);
    error(:,:,i) = readmatrix([path, MethodName, '\errorX.csv']);
    i = i + 1;
    
    DataIndex = 32; %RND

    DataToPlot = data(:,DataIndex, mi); %run/step
    Temprature = data(:,1, mi);
    % h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(Temprature, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 4
         h.MarkerFaceColor = h.Color;
    % end
    hold on
    legendStrings = [legendStrings, string([ MethodName])];;

end

% xlim([1e-3 1e1])
% ylim([1e1 2e5])
yl = ylim;
xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
xticks([1e-3, 1e-2, 1e-1, 1e0])
yticks([1e2, 1e3, 1e4])


myFontSize =12;
set(gca,'fontsize',myFontSize) 
ylabel ('$N_\mathrm{RND,acc}$', 'FontSize', myFontSize)
xlabel ('$T$', 'FontSize', myFontSize)
% leg1 = legend(legendStrings, 'Location','best','Interpreter','latex', 'FontSize', myFontSize);
% set(leg1,'Box','off')

set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
% pbaspect([2 1 1])
set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot

% exportgraphics(gcf,[save_address '/', 'DythinVsFukui_RND_accept.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'DythinVsFukui_RND_accept.eps'])

%% save
% set(axs(1),'XTickLabel','');
% % set(axs(2:2:end),'YTickLabel','')
% xlabel(axs(1), '');

linkaxes([axs(1:2)],'x')
ListSubfigLabels = {'(a)', '(b)', '(c)'}
AddLetters2Plots(fig, ListSubfigLabels, 'HShift', 0.02, 'VShift', 0, ...
    'Direction', 'TopDown', 'FontWeight', 'normal', 'FontSize', myFontSize)

filename = ['DythinVsFukui_Temperature'];
exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])

%% load data Sizes
clc
clear all
path = ['E:\visualize\NewAnalysis\DynamicThinning_Vs_Fukui\Sizes\'];
methods = ["DyThin", "FT"];


% MethodName = sprintf('%s',method);
% data(:,:,i) =  readmatrix([path, MethodName, '\MeanX.csv']);
% error(:,:,i) = readmatrix([path, MethodName, '\errorX.csv']);


date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
save_address = ['E:/visualize/output/', date_time_string , '_DythinVsFukui'];
[status, msg, msgID] = mkdir(save_address)
dpi = 600;

listL = [16, 24, 32, 48, 64];
%%---------------------------------------------------------------------------------
%% plot C sizes

f = figure;
f.Position = [100 100 300 200];
cmap = linspecer(length(methods));
markerShapes = {'v', 's'};
legendStrings=[];
mi=0;
for method = methods
    data=[];
    error=[];
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    i = 1;
    for L=listL
        Lstr = num2str(L);
        data_path = [path, MethodName, '\' , Lstr , '\'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
        i = i + 1;
    end
    
    DataIndex = 31; %Pair

    DataToPlot=data(DataIndex,:);
    % h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 4
         h.MarkerFaceColor = h.Color;
    % end
    hold on
    legendStrings = [legendStrings, string([ MethodName])];;

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
ylabel ('Pairwise Energy/Spin', 'FontSize', myFontSize)
xlabel ('$L$', 'FontSize', myFontSize)
leg1 = legend(legendStrings, 'Location','east','Interpreter','latex', 'FontSize', myFontSize);
set(leg1,'Box','off')

set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
% pbaspect([2 1 1])

exportgraphics(gcf,[save_address '/', 'DythinVsFukui_C_Tc.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/', 'DythinVsFukui_C_Tc.eps'])

%%---------------------------------------------------------------------------------
%% plot RND sizes

f = figure;
f.Position = [100 100 300 200];
cmap = linspecer(length(methods));
markerShapes = {'v', 's'};
legendStrings=[];

mi=0;
for method = methods
    data=[];
    error=[];
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    i = 1;
    for L=listL
        Lstr = num2str(L);
        data_path = [path, MethodName, '\' , Lstr , '\'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
        i = i + 1;
    end
    
    DataIndex = 32; %RND

    DataToPlot=data(DataIndex,:);
    % h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 4
         h.MarkerFaceColor = h.Color;
    % end
    hold on
    legendStrings = [legendStrings, string([ MethodName])];;

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
ylabel ('Random Number/Spin', 'FontSize', myFontSize)
xlabel ('$L$', 'FontSize', myFontSize)
leg1 = legend(legendStrings, 'Location','west','Interpreter','latex', 'FontSize', myFontSize);
set(leg1,'Box','off')

set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
% pbaspect([2 1 1])

exportgraphics(gcf,[save_address '/', 'DythinVsFukui_RND_Tc.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/', 'DythinVsFukui_RND_Tc.eps'])


