%% load data
clc
clear all


date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
save_address = ['E:/visualize/output/', date_time_string , '_Tomita'];
[status, msg, msgID] = mkdir(save_address)
dpi = 600;
%% t_eff Plain Main

f = figure;
f.Position = [100 100 300 200];

markerShapes = {'^', 'v', 's', 'o', 'v'};

path = ['E:\visualize\NewAnalysis\Tomita different AlphaTilde\'];
listAlphaString = readmatrix([path, 'list.txt'],'OutputType' , 'string');
listAlpha = readmatrix([path, 'list.txt']);
N = length(listAlpha);
data =[];
error=[];
mi = 1;
for i = 1:length(listAlphaString)
    AlphaStr = listAlphaString(i);
    Alphachar = sprintf('%s',AlphaStr);
    data_path = [path, Alphachar , '/'];

    % [status, msg, msgID] = mkdir(folder_address)

    data(:,i) = readmatrix([data_path 'MeanX.csv']);
    error(:,i) = readmatrix([data_path 'errorX.csv']);
end

%Pair
DataIndex = 31; filename='Pair'; Ylab='$C_{\mathrm{pairwise}}$'; %% per spin
% RND
% DataIndex = 32; filename='RND'; Ylab='Random Number/Spin'; %% per spin
%t_MCS
% DataIndex = 30; filename='t_MCS'; Ylab='$t_{\mathrm{MCS}}$'; %% per sweep
%Accept
% DataIndex = 26; filename='Accept'; Ylab='Acceptance Ratio';


% h = errorbar(listAlpha, data(DataIndex,:), error(DataIndex,:), 'Marker', 'o' ,'MarkerSize', 7, ...
%         'LineWidth', 0.8);
h = plot(listAlpha, data(DataIndex,:), 'Marker', markerShapes{mi}, 'Marker', 'o' ,'MarkerSize', 5, ...
        'LineWidth', 0.5);


h.MarkerFaceColor = h.Color;
% end
% set(h, {'MarkerFaceColor'}, get(h,'Color')); 
hold on


% ylim([1e2 1e6])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =12;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
ylabel (Ylab)

xlabel('$\tilde{\alpha}$', 'FontSize', myFontSize)
% leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', myFontSize);
% set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot

% pbaspect([2 1 1])
% 

filename = ['TomitaAlphaTilde', filename];
exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])


%%---------------------------------------------------------------------------------
%% t_MCS Plain Main

clf('reset')


methods=["clock", "Tomita"];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'^', 's'};
legendStrings=[];
tauExponents=[];
mi=0;
for method = methods
    mi = mi+1;
    
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = ['E:\visualize\NewAnalysis\Tomita AlphaTilde=0/' MethodName '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = num2str(listL(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    
    %Pair
    % DataIndex = 31; filename='Pair'; Ylab='Pairwise Energy/Spin'; %% per spin
    % RND
    % DataIndex = 32; filename='RND'; Ylab='Random Number/Spin'; %% per spin
    %t_MCS
    % DataIndex = 30; filename='t_MCS'; Ylab='$t_{\mathrm{MCS}}$'; %% per sweep
    %Accept
    DataIndex = 26; filename='Accept'; Ylab='Acceptance Ratio';

    DataToPlot=data(DataIndex,:); %run/step
    ErrorOfPlot = error(DataIndex,:);
    
    h = errorbar(listL, DataToPlot, ErrorOfPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    % h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);

    % if mi >= 4
    h.MarkerFaceColor = h.Color;
    % end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
    if strcmp(method, "Tomita")
        MethodName = 'Tomita ($\tilde{\alpha}=0$)';
    end
    legendStrings = [legendStrings, string([ MethodName])];;

end

% ylim([1e2 1e6])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =12;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
ylabel (Ylab)
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$', 'FontSize', myFontSize)

xlabel('$L$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
leg1 = legend(legendStrings, 'Location','southwest','Interpreter','latex', 'FontSize', myFontSize);
set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot

pbaspect([2 1 1])
% 
filename = ['Tomita alphaTilde=0 Vs Clock ', filename];
exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])
