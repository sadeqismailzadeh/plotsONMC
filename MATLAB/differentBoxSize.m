%% load data
clc
clear all
path = ['E:\visualize\NewAnalysis\Different Box\'];
methods = ["clock", "SCO", "Tomita"];


date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
save_address = ['E:/visualize/output/', date_time_string , '_DifferentBoxing'];
[status, msg, msgID] = mkdir(save_address)

listboxes = {'0' '0.550' '0.650' '0.700' '0.750' ...
    '0.800' '0.850' '0.900' '0.950' '0.970' '0.980' '0.990'};
listboxesNumeric = [0 0.55  0.65 0.7 0.75 0.8 0.85 0.9 0.95 ...
    0.97 0.98 0.99];
dpi = 600;
%%---------------------------------------------------------------------------------
%% plot C

clf('reset')
cmap = linspecer(length(methods));
markerShapes = {'v', 's', '^'};
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
    for Box=listboxes
        % Boxstr = num2str(Box);
        Boxstr = Box{1};
        data_path = [path, MethodName, '\' , Boxstr , '\'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
        i = i + 1;
    end

    %Pair
    % DataIndex = 31; filename='Pair'; Ylab='Pairwise Energy/Spin'; %% per spin
    %RND
    % DataIndex = 32; filename='RND'; Ylab='Random Number/Spin'; %% per spin
    %t_MCS
    % DataIndex = 30; filename='t_sweep'; Ylab='$t_{\mathrm{MCS}}$'; %% per sweep
    %Accept
    DataIndex = 26; filename='Accept'; Ylab='Acceptance Ratio';
    
    DataToPlot=data(DataIndex,:);
    % h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    % h = plot(listboxesNumeric, DataToPlot);
    h = plot(listboxesNumeric, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    % if strcmp(method, "Over-relaxed Clock")
    %     h.LineWidth = 2.5;
    %     h.MarkerSize = 9;
    % end
    h.MarkerFaceColor = h.Color;

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
ylabel (Ylab, 'FontSize', myFontSize)
xlabel ('$E_{\mathrm{group}} / E_\infty$', 'FontSize', myFontSize)
% leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', myFontSize);
% leg1 = legend(legendStrings, 'Location','west','Interpreter','latex', 'FontSize', myFontSize);
% leg1 = legend(legendStrings, 'Location','northeast','Interpreter','latex', 'FontSize', myFontSize);
leg1 = legend(legendStrings, 'Location','southeast','Interpreter','latex', 'FontSize', myFontSize);

set(leg1,'Box','off')

set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

% set(gca, 'YScale', 'log')   %using plot
pbaspect([2 1 1])
filename = ['boxEproportion',filename];
exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])