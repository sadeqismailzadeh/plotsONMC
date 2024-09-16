%% load data
clc
clear all
% path = ['E:\visualize\NewAnalysis\Different Box vs L\'];
% path = ['E:\visualize\NewAnalysis\SCO\Validity\Overrelaxed\'];


date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
save_address = ['E:/visualize/output/', date_time_string , '_Box_r_Vs_L'];
[status, msg, msgID] = mkdir(save_address)
dpi = 600;

%%---------------------------------------------------------------------------------
%% new subfig
fig = figure;
fig.Position  = [0 0 450 600];
% tiledlayout(4,2)
axs = tight_subplot(4,2,[.03 .05],[.09 .03],[.15 .03])

%%---------------------------------------------------------------------------------
%% t_MCS Plain Main subfigs

% clf('reset')
for ii=1:8
% methods=["clock", "SCO" , "SCO (OR)" , "Tomita"];
methods=["32", "64", "128" , "256"];
SystemLength=[32, 64, 128 , 256];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];
subpathclock = 'E:\visualize\NewAnalysis\Different Box vs L i7 6700K\random\clock';
subpathSCO = 'E:\visualize\NewAnalysis\Different Box vs L i7 6700K\random\SCO';
subpathTomita = 'E:\visualize\NewAnalysis\Different Box vs L i7 6700K\random\Tomita';

if (mod(ii,2) == 0)
    subpath = subpathSCO;
else
    subpath = subpathclock;
end

cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', 'o', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};
legendStrings=[];
tauExponents=[];
mi=0;
axes(axs(ii))
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = [subpath '\' MethodName '/'];
    listr = readmatrix([path, 'list.txt']);
    listrStr = readmatrix([path, 'list.txt'], 'OutputType', 'string');
    listN = listr.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = sprintf('%s',listrStr(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end

    Nsystem = 1;
    if (ii <= 2)
        %Pair
        DataIndex = 31; filename='Pair'; Ylab='$N_\mathrm{PW}$'; %% per spin
    elseif (ii > 2 && ii <=4)
        % RND
        DataIndex = 32; filename='RND'; Ylab='$N_\mathrm{RND}$'; %% per spin
    elseif (ii > 4 && ii <=6)
        % t_MCS
        DataIndex = 30; filename='t_MCS'; Ylab='$t_{\mathrm{spin}} \; [\mu s]$'; Nsystem = SystemLength(mi) ^2 * 1e-6; %% per
        % spin
    elseif (ii > 6 && ii <=8)
        %Accept
        DataIndex = 26; filename='Accept'; Ylab='$P_{\mathrm{acc}}$';
    end

    %Accept OR
    % DataIndex = 27; filename='Accept_OR'; Ylab='OR Acceptance Ratio';

    DataToPlot= data(DataIndex,:) / Nsystem; %run/step
    ErrorOfPlot=error(DataIndex,:) / Nsystem;
    % h = errorbar(listr, DataToPlot, ErrorOfPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    
    h = plot(listr, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 5, ...
            'Color',cmap(mi,:), 'LineWidth', 0.6);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 4
    h.MarkerFaceColor = h.Color;
    % end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    legendStrings = [legendStrings, string(['$L = ' MethodName, '$'])];

end

% ylim([2e-7 10e-7])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =11;
set(gca,'fontsize',myFontSize) 




% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','northeast','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','north','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','southeast','Interpreter','latex', 'FontSize', 10);

if (ii == 1)
    leg1 = legend(legendStrings, 'Location','best','Interpreter','latex', 'FontSize', 10);
    set(leg1,'Box','off')
end


if (ii <= 2)
    set(gca, 'YScale', 'log')   %using plot
end

if (ii == 1)
    title ('clock')
end

if (ii == 2)
    title ('SCO')
end


% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

ylabel (Ylab, 'FontSize', myFontSize)   % show y label
xlabel('$r_\mathrm{box}$', 'FontSize', myFontSize) % show x label



% set(ax,'xticklabel',[])  % hide x axis  % show y label
% set(ax,'yticklabel',[])  % hide y axis  % show x label


% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot

% pbaspect([1.5 1 1])
% 
end
%%---------------------------------------------------------------------------------
%% remove labels names
set(axs(1:end-2),'XTickLabel','');
set(axs(2:2:end),'YTickLabel','')
xlabel(axs(1:end-2), '');
ylabel(axs(2:2:end), '');

for jj=1:4
beginIndex = 2*jj-1;
endIndex = 2*jj;
linkaxes([axs(beginIndex:endIndex)],'xy')
end
% filename = ['Box_r_Vs_L_clock_', filename];



%%---------------------------------------------------------------------------------
%% add labels to subfigs
ListSubfigLabels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)'}
AddLetters2Plots(fig, ListSubfigLabels, 'HShift', 0.31, 'VShift', 0.028, ...
    'Direction', 'TopDown', 'FontWeight', 'normal', 'FitLocation', true)
%%---------------------------------------------------------------------------------
%% save

% filename = ['Box_r_Vs_L_clock_OR_', filename];
% filename = ['Box_r_Vs_L_clock_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_', filename];
% filename = ['Box_r_Vs_L_SCO_OR_', filename];

filename = 'Box_r_Vs_L';
exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])




%%---------------------------------------------------------------------------------
%% new subfig L=256
fig = figure;
fig.Position  = [0 0 500 330];
% tiledlayout(4,2)
axs = tight_subplot(2,2,[.05 .13],[.15 .03],[.12 .03])

%%---------------------------------------------------------------------------------
%% L=256 t_MCS Plain Main subfig

% clf('reset')
for ii=1:4
% methods=["clock", "SCO" , "SCO (OR)" , "Tomita"];
% methods=["32", "64", "128" , "256"];
% SystemLength=[32, 64, 128 , 256];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];
subpathclock = 'E:\visualize\NewAnalysis\Different Box vs L i7 6700K\random\clock';
subpathSCO = 'E:\visualize\NewAnalysis\Different Box vs L i7 6700K\random\SCO';
subpathTomita = 'E:\visualize\NewAnalysis\Different Box vs L i7 6700K\random\Tomita';

subpathclockOR = 'E:\visualize\NewAnalysis\Different Box vs L i7 6700K\OR\clock';
subpathSCOOR = 'E:\visualize\NewAnalysis\Different Box vs L i7 6700K\OR\SCO';
subpathTomitaOR = 'E:\visualize\NewAnalysis\Different Box vs L i7 6700K\OR\Tomita';

methods=["clock", "SCO" , "Tomita"];

isOR = false;
if (isOR)
    subpaths={subpathclockOR, subpathSCOOR, subpathTomitaOR}; % overrelaxed
else
    subpaths={subpathclock, subpathSCO, subpathTomita};  %random
end



cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', 'o', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};
legendStrings=[];
tauExponents=[];
mi=0;
axes(axs(ii))
SystemLength=256;
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = [subpaths{mi}, '/' num2str(SystemLength) '/'];  % for reference version
    listr = readmatrix([path, 'list.txt']);
    listrStr = readmatrix([path, 'list.txt'], 'OutputType', 'string');
    listN = listr.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = sprintf('%s',listrStr(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end

    Nsystem = 1;
    if (ii ==1)
        %Pair
        DataIndex = 31; filename='Pair'; Ylab='$N_\mathrm{pw}$'; %% per spin
    elseif (ii == 3)
        % RND
        DataIndex = 32; filename='RND'; Ylab='$N_\mathrm{RND}$'; %% per spin
    elseif (ii == 2)
        % t_MCS
        DataIndex = 30; filename='t_MCS'; Ylab='$t_{\mathrm{spin}} \; [\mu s]$'; Nsystem = SystemLength ^2 * 1e-6; %% per
        % spin
    elseif (ii == 4)
        if (isOR)
            %Accept OR
            DataIndex = 27; filename='Accept_OR'; Ylab='$P_{\mathrm{acc, OR}}$';
        else
            %Accept
            DataIndex = 26; filename='Accept'; Ylab='$P_{\mathrm{acc}}$';
        end
        
    end

    

    DataToPlot= data(DataIndex,:) / Nsystem; %run/step
    ErrorOfPlot=error(DataIndex,:) / Nsystem;
    % h = errorbar(listr, DataToPlot, ErrorOfPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listr, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 6, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 4
    h.MarkerFaceColor = h.Color;
    % end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    legendStrings = [legendStrings, string([ MethodName])];
end

% ylim([2e-7 10e-7])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =13;
set(gca,'fontsize',myFontSize) 




% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','northeast','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','north','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','southeast','Interpreter','latex', 'FontSize', 10);

if (ii == 1)
    leg1 = legend(legendStrings, 'Location','southeast','Interpreter','latex', 'FontSize', myFontSize);
    set(leg1,'Box','off')
    set(gca, 'YScale', 'log')   %using plot
elseif (ii == 2)
    if (isOR)
        %Accept OR
    else
        %Accept
    end
elseif (ii == 3)
    set(gca, 'YScale', 'log')   %using plot
    if (isOR)
        %Accept OR
    else
        yticks([1e1, 1e2])
    end
elseif (ii == 4)
    if (isOR)
        yticks([ 0.3, 0.6, 0.9])
    else
        %Accept
    end
    
end


% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

ylabel (Ylab, 'FontSize', myFontSize)   % show y label
xlabel('$r_\mathrm{box}$', 'FontSize', myFontSize) % show x label



% set(ax,'xticklabel',[])  % hide x axis  % show y label
% set(ax,'yticklabel',[])  % hide y axis  % show x label


% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot

% pbaspect([1.5 1 1])
% 
end
%%---------------------------------------------------------------------------------
%% remove labels names
set(axs(1:end-2),'XTickLabel','');
% % set(axs(2:2:end),'YTickLabel','')
xlabel(axs(1:end-2), '');
% ylabel(axs(2:2:end), '');
linkaxes([axs(1:end)],'x')
% for jj=1:4
% beginIndex = 2*jj-1;
% endIndex = 2*jj;
% linkaxes([axs(beginIndex:endIndex)],'xy')
% end
% filename = ['Box_r_Vs_L_clock_', filename];



%%---------------------------------------------------------------------------------
%% add labels to subfigs
ListSubfigLabels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)'}

if (isOR)
AddLetters2Plots(fig, ListSubfigLabels, 'HShift', 0.02, 'VShift', 0.0, ...
    'Direction', 'TopDown', 'FontWeight', 'normal', 'FitLocation', true, 'FontSize', myFontSize)
else
AddLetters2Plots(fig, ListSubfigLabels, 'HShift', 0.28, 'VShift', 0.05, ...
    'Direction', 'TopDown', 'FontWeight', 'normal', 'FitLocation', true, 'FontSize', myFontSize)
end


%%---------------------------------------------------------------------------------
%% save

% filename = ['Box_r_Vs_L_clock_OR_', filename];
% filename = ['Box_r_Vs_L_clock_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_', filename];
% filename = ['Box_r_Vs_L_SCO_OR_', filename];
if (isOR)
    filename = 'Box_r_Vs_L_compare_OR_L256';
else
    filename = 'Box_r_Vs_L_compare_L256';
end

exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])




%%---------------------------------------------------------------------------------
%% new subfig clock multiple sizes
fig = figure;
fig.Position  = [0 0 500 330];
% tiledlayout(4,2)
axs = tight_subplot(2,2,[.05 .13],[.15 .03],[.12 .03])

%%---------------------------------------------------------------------------------
%% clock multiple sizes Plain Main subfig

% clf('reset')
for ii=1:4
% methods=["clock", "SCO" , "SCO (OR)" , "Tomita"];
methods=["32", "64", "128" , "256"];
SystemLength=[32, 64, 128 , 256];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];
subpathclock = 'E:\visualize\NewAnalysis\Different Box vs L i7 6700K\random\clock';
subpathSCO = 'E:\visualize\NewAnalysis\Different Box vs L i7 6700K\random\SCO';
subpathTomita = 'E:\visualize\NewAnalysis\Different Box vs L i7 6700K\random\Tomita';

subpath = subpathclock;

isOR=false;


cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', 'o', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};
legendStrings=[];
tauExponents=[];
mi=0;
axes(axs(ii))
for method = methods
     mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = [subpath '\' MethodName '/'];
    listr = readmatrix([path, 'list.txt']);
    listrStr = readmatrix([path, 'list.txt'], 'OutputType', 'string');
    listN = listr.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = sprintf('%s',listrStr(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end

    Nsystem = 1;
    if (ii ==1)
        %Pair
        DataIndex = 31; filename='Pair'; Ylab='$N_\mathrm{pw}$'; %% per spin
    elseif (ii == 3)
        % RND
        DataIndex = 32; filename='RND'; Ylab='$N_\mathrm{RND}$'; %% per spin
    elseif (ii == 2)
        % t_MCS
        DataIndex = 30; filename='t_MCS'; Ylab='$t_{\mathrm{spin}} \; [\mu s]$'; Nsystem = SystemLength(mi) ^2 * 1e-6; %% per
        % spin
    elseif (ii == 4)
        if (isOR)
            %Accept OR
            DataIndex = 27; filename='Accept_OR'; Ylab='$P_{\mathrm{acc, OR}}$';
        else
            %Accept
            DataIndex = 26; filename='Accept'; Ylab='$P_{\mathrm{acc}}$';
        end
        
    end

    

    DataToPlot= data(DataIndex,:) / Nsystem; %run/step
    ErrorOfPlot=error(DataIndex,:) / Nsystem;
    % h = errorbar(listr, DataToPlot, ErrorOfPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listr, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 6, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 4
    h.MarkerFaceColor = h.Color;
    % end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    legendStrings = [legendStrings, string(['$L = ' MethodName, '$'])];
end

% ylim([2e-7 10e-7])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =13;
set(gca,'fontsize',myFontSize) 




% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','northeast','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','north','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','southeast','Interpreter','latex', 'FontSize', 10);



if (ii == 1)
    leg1 = legend(legendStrings, 'Location','southeast','Interpreter','latex', 'FontSize', myFontSize-1);
    set(leg1,'Box','off')
    set(gca, 'YScale', 'log')   %using plot
elseif (ii == 2)
    %
elseif (ii == 3)
    yticks([4, 8, 12])
    set(gca, 'YScale', 'log')   %using plot
elseif (ii == 4)
    
end

% if (ii == 1)
%     title ('clock')
% end
% 
% if (ii == 2)
%     title ('SCO')
% end


% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

ylabel (Ylab, 'FontSize', myFontSize)   % show y label
xlabel('$r_\mathrm{box}$', 'FontSize', myFontSize) % show x label



% set(ax,'xticklabel',[])  % hide x axis  % show y label
% set(ax,'yticklabel',[])  % hide y axis  % show x label


% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot

% pbaspect([1.5 1 1])
% 
end
%%---------------------------------------------------------------------------------
%% remove labels names
set(axs(1:end-2),'XTickLabel','');
% % set(axs(2:2:end),'YTickLabel','')
xlabel(axs(1:end-2), '');
% ylabel(axs(2:2:end), '');
linkaxes([axs(1:end)],'x')
% for jj=1:4
% beginIndex = 2*jj-1;
% endIndex = 2*jj;
% linkaxes([axs(beginIndex:endIndex)],'xy')
% end
% filename = ['Box_r_Vs_L_clock_', filename];



%%---------------------------------------------------------------------------------
%% add labels to subfigs
ListSubfigLabels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)'}


AddLetters2Plots(fig, ListSubfigLabels, 'HShift', 0.015, 'VShift', 0.0, ...
    'Direction', 'TopDown', 'FontWeight', 'normal', 'FitLocation', true, 'FontSize', myFontSize)

% AddLetters2Plots(fig, ListSubfigLabels, 'HShift', 0.28, 'VShift', 0.05, ...
%     'Direction', 'TopDown', 'FontWeight', 'normal', 'FitLocation', true)



%%---------------------------------------------------------------------------------
%% save

% filename = ['Box_r_Vs_L_clock_OR_', filename];
% filename = ['Box_r_Vs_L_clock_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_', filename];
% filename = ['Box_r_Vs_L_SCO_OR_', filename];

filename = 'Box_r_Vs_L_clock';


exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])


%%---------------------------------------------------------------------------------
%% t_MCS Plain Main Single

f = figure;
f.Position = [100 100 300 200];
% clf('reset')
% methods=["clock", "SCO" , "SCO (OR)" , "Tomita"];
methods=["32", "64", "128" , "256"];
SystemLength=[32, 64, 128 , 256];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];
subpath1 = 'E:\visualize\NewAnalysis\Different Box vs L\SCO (equilibrated)';
subpath2 = 'E:\visualize\NewAnalysis\Different Box vs L\clock square (uneqilibrated)';
subpath3 = 'E:\visualize\NewAnalysis\Different Box vs L\clock (equilibrated)';
subpath4 = 'E:\visualize\NewAnalysis\Different Box vs L\clock DyThin (equilibrated)';
subpath5 = 'E:\visualize\NewAnalysis\Different Box vs L\SCO DyThin (equilibrated)';
subpath6 = 'E:\visualize\NewAnalysis\Different Box vs L\Tomita (equilibrated)';
subpath7 = 'E:\visualize\NewAnalysis\Different Box vs L\clock OR (equilibrated)';
subpath8 = 'E:\visualize\NewAnalysis\Different Box vs L\SCO OR (equilibrated)';
subpath9 = 'E:\visualize\NewAnalysis\Different Box vs L\clock (unequilibrated)';
subpath10 = 'E:\visualize\NewAnalysis\Different Box vs L\SCO (equilibrated) 2';
subpath11 = 'E:\visualize\NewAnalysis\Different Box vs L\Tomita OR (equilibrated)';

subpath = subpath3;
cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', 'o', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};
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
    path = [subpath '\' MethodName '/'];
    listr = readmatrix([path, 'list.txt']);
    listrStr = readmatrix([path, 'list.txt'], 'OutputType', 'string');
    listN = listr.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = sprintf('%s',listrStr(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end

    Nsystem = 1;
    %Pair
    % DataIndex = 31; filename='Pair'; Ylab='$N_\mathrm{pairwise} / \mathrm{spin}$'; %% per spin
    % % RND
    % DataIndex = 32; filename='RND'; Ylab='$N_\mathrm{RND}/ \mathrm{spin}$'; %% per spin
    % t_MCS
    % DataIndex = 30; filename='t_MCS'; Ylab='$t_{\mathrm{MCS}} / \mathrm{spin}$'; Nsystem = SystemLength(mi) ^2; %% per
    % spin
    %Accept
    DataIndex = 26; filename='Accept'; Ylab='$P_{\mathrm{acc}}$';
    %Accept OR
    % DataIndex = 27; filename='Accept_OR'; Ylab='OR Acceptance Ratio';

    DataToPlot= data(DataIndex,:) / Nsystem; %run/step
    ErrorOfPlot=error(DataIndex,:) / Nsystem;
    % h = errorbar(listr, DataToPlot, ErrorOfPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    
    h = plot(listr, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 5, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 4
    h.MarkerFaceColor = h.Color;
    % end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    legendStrings = [legendStrings, string(['$L = ' MethodName, '$'])];

end

% ylim([2e-7 10e-7])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =12;
set(gca,'fontsize',myFontSize) 




% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','northeast','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','north','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','southeast','Interpreter','latex', 'FontSize', 10);
leg1 = legend(legendStrings, 'Location','best','Interpreter','latex', 'FontSize', 10);

set(leg1,'Box','off')


% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

ylabel (Ylab, 'FontSize', myFontSize)   % show y label
xlabel('$r_\mathrm{box}$', 'FontSize', myFontSize) % show x label



% set(ax,'xticklabel',[])  % hide x axis  % show y label
% set(ax,'yticklabel',[])  % hide y axis  % show x label


% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot

% pbaspect([1.5 1 1])
filename = ['Box_r_Vs_L_clock_', filename];
% filename = ['Box_r_Vs_L_clock_OR_', filename];
% filename = ['Box_r_Vs_L_clock_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_', filename];
% filename = ['Box_r_Vs_L_SCO_OR_', filename];
exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])


%%---------------------------------------------------------------------------------
%% Delta t_MCS Plain Main

% clf('reset')
f = figure;
f.Position = [100 100 300 200];
% methods=["clock", "SCO" , "SCO (OR)" , "Tomita"];
methods=["32", "64", "128" , "256"];
SystemLength=[32, 64, 128 , 256];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];
subpath1 = 'E:\visualize\NewAnalysis\Different Box vs L\SCO (equilibrated)';
subpath2 = 'E:\visualize\NewAnalysis\Different Box vs L\clock square (uneqilibrated)';
subpath3 = 'E:\visualize\NewAnalysis\Different Box vs L\clock (equilibrated)';
subpath4 = 'E:\visualize\NewAnalysis\Different Box vs L\clock DyThin (equilibrated)';
subpath5 = 'E:\visualize\NewAnalysis\Different Box vs L\SCO DyThin (equilibrated)';
subpath6 = 'E:\visualize\NewAnalysis\Different Box vs L\Tomita (equilibrated)';
subpath7 = 'E:\visualize\NewAnalysis\Different Box vs L\clock OR (equilibrated)';
subpath8 = 'E:\visualize\NewAnalysis\Different Box vs L\SCO OR (equilibrated)';
subpath9 = 'E:\visualize\NewAnalysis\Different Box vs L\clock (unequilibrated)';
subpath10 = 'E:\visualize\NewAnalysis\Different Box vs L\SCO (equilibrated) 2';
subpath11 = 'E:\visualize\NewAnalysis\Different Box vs L\Tomita OR (equilibrated)';

subpath = subpath9;
cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', '^', 'v'};
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
    path = [subpath10 '\' MethodName '/'];  % for reference version
    listr = readmatrix([path, 'list.txt']);
    listrStr = readmatrix([path, 'list.txt'], 'OutputType', 'string');
    listN = listr.^2;
    N = length(listN)
    data1 =[];
    error1=[];
    for i=1:N
        name = sprintf('%s',listrStr(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data1(:,i) = readmatrix([data_path 'MeanX.csv']);
        error1(:,i) = readmatrix([data_path 'errorX.csv']);
    end

  
    Nsystem = 1;
    %Pair
    % DataIndex = 31; filename='Pair'; Ylab='$C_\mathrm{pairwise}$'; %% per spin
    % % RND
    % DataIndex = 32; filename='RND'; Ylab='$C_\mathrm{RND}$'; %% per spin
    % t_MCS
    DataIndex = 30; filename='t_MCS'; Ylab='$t_{\mathrm{update}}$'; Nsystem = SystemLength(mi) ^2; %% per
    % spin
    %Accept
    % DataIndex = 26; filename='Accept'; Ylab='Acceptance Ratio';
    %Accept OR
    % DataIndex = 27; filename='Accept_OR'; Ylab='OR Acceptance Ratio';

    DataToPlot= data(DataIndex,:) / Nsystem; %run/step
    ErrorOfPlot=error(DataIndex,:) / Nsystem;
    % h = errorbar(listr, DataToPlot, ErrorOfPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listr, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 5, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 4
    h.MarkerFaceColor = h.Color;
    % end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    legendStrings = [legendStrings, string(['$L = ' MethodName, '$'])];

end
xlim([0, 35])
% ylim([2e-7 10e-7])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =12;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
ylabel (Ylab, 'FontSize', myFontSize)
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$', 'FontSize', myFontSize)

xlabel('$r_\mathrm{box}$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','northeast','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','north','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','southeast','Interpreter','latex', 'FontSize', 10);
leg1 = legend(legendStrings, 'Location','best','Interpreter','latex', 'FontSize', 10);

set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
% axis 'padded'
% pbaspect([1.5 1 1])
% 
% filename = ['Box_r_Vs_L_clock_', filename];
% filename = ['Box_r_Vs_L_clock_OR_', filename];
% filename = ['Box_r_Vs_L_clock_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_', filename];
% filename = ['Box_r_Vs_L_Diff_fukui_Dythin_', filename];
filename = ['Box_r_Vs_L_compare_', filename];
exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])


%%---------------------------------------------------------------------------------
%% L=256 t_MCS Plain Main

% clf('reset')
f = figure;
f.Position = [100 100 300 200];
% methods=["clock", "SCO" , "SCO (OR)" , "Tomita"];
% methods=["32", "64", "128" , "256"];
% SystemLength=[32, 64, 128 , 256];
% methods=["null", "null", "null", "OR-Clock", "null"];

% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];
subpath1 = 'E:\visualize\NewAnalysis\Different Box vs L\SCO (equilibrated)';
subpath2 = 'E:\visualize\NewAnalysis\Different Box vs L\clock square (uneqilibrated)';
subpath3 = 'E:\visualize\NewAnalysis\Different Box vs L\clock (equilibrated)';
subpath4 = 'E:\visualize\NewAnalysis\Different Box vs L\clock DyThin (equilibrated)';
subpath5 = 'E:\visualize\NewAnalysis\Different Box vs L\SCO DyThin (equilibrated)';
subpath6 = 'E:\visualize\NewAnalysis\Different Box vs L\Tomita (equilibrated)';
subpath7 = 'E:\visualize\NewAnalysis\Different Box vs L\clock OR (equilibrated)';
subpath8 = 'E:\visualize\NewAnalysis\Different Box vs L\SCO OR (equilibrated)';
subpath9 = 'E:\visualize\NewAnalysis\Different Box vs L\clock (unequilibrated)';
subpath10 = 'E:\visualize\NewAnalysis\Different Box vs L\SCO (equilibrated) 2';
subpath11 = 'E:\visualize\NewAnalysis\Different Box vs L\Tomita OR (equilibrated)';

% methods=["clock", "SCO" , "Tomita"];
methods=["clock", "clock non-eq"];
subpaths={subpath3, subpath9}; % non-eq

% subpaths={subpath3, subpath10, subpath6}; % Plain
% subpaths={subpath7, subpath8, subpath11}; % overrelaxed

cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', '^', 'v'};
legendStrings=[];
tauExponents=[];
mi=0;
SystemLength=256;
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = [subpaths{mi}, '/' num2str(SystemLength) '/'];  % for reference version
    listr = readmatrix([path, 'list.txt']);
    listrStr = readmatrix([path, 'list.txt'], 'OutputType', 'string');
    listN = listr.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = sprintf('%s',listrStr(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end

    Nsystem = 1;
    %Pair
    % DataIndex = 31; filename='Pair'; Ylab='$C_\mathrm{pairwise}$'; %% per spin
    % % RND
    % DataIndex = 32; filename='RND'; Ylab='$C_\mathrm{RND}$'; %% per spin
    % t_MCS
    % DataIndex = 30; filename='t_MCS'; Ylab='$t_{\mathrm{update}}$'; Nsystem = SystemLength ^2; %% per
    % spin
    %Accept
    DataIndex = 26; filename='Accept'; Ylab='Acceptance Ratio';
    %Accept OR
    % DataIndex = 27; filename='Accept_OR'; Ylab='OR Acceptance Ratio';

    DataToPlot= data(DataIndex,:) / Nsystem; %run/step
    ErrorOfPlot=error(DataIndex,:) / Nsystem;
    % h = errorbar(listr, DataToPlot, ErrorOfPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listr, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 5, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 4
    h.MarkerFaceColor = h.Color;
    % end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    legendStrings = [legendStrings, string([ MethodName])];

end
xlim([0, 35])
% ylim([2e-7 10e-7])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =12;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
ylabel (Ylab, 'FontSize', myFontSize)
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$', 'FontSize', myFontSize)

xlabel('$r_\mathrm{box}$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','northeast','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','north','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','southeast','Interpreter','latex', 'FontSize', 10);
leg1 = legend(legendStrings, 'Location','best','Interpreter','latex', 'FontSize', 10);

set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
% axis 'padded'
% pbaspect([1.5 1 1])
% 
% filename = ['Box_r_Vs_L_clock_', filename];
% filename = ['Box_r_Vs_L_clock_OR_', filename];
% filename = ['Box_r_Vs_L_clock_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_', filename];
% filename = ['Box_r_Vs_L_Diff_fukui_Dythin_', filename];
% filename = ['Box_r_Vs_L_compare_', filename];
filename = ['Box_r_Vs_L_compare_OR_L256_', filename];
exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])



%%---------------------------------------------------------------------------------
%% L=256 OR=10 t_MCS Main

% clf('reset')
f = figure;
f.Position = [100 100 300 200];
% methods=["clock", "SCO" , "SCO (OR)" , "Tomita"];
% methods=["32", "64", "128" , "256"];
% SystemLength=[32, 64, 128 , 256];
% methods=["null", "null", "null", "OR-Clock", "null"];

% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];
subpath1 = 'E:\visualize\NewAnalysis\r for minimum t_MCS OR=10 L=256\clock';
subpath2 = 'E:\visualize\NewAnalysis\r for minimum t_MCS OR=10 L=256\SCO';
subpath3 = 'E:\visualize\NewAnalysis\r for minimum t_MCS OR=10 L=256\Tomita';

methods=["clock", "SCO" , "Tomita"];

subpaths={subpath1, subpath2, subpath3};

cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', '^', 'v'};
legendStrings=[];
tauExponents=[];
mi=0;
SystemLength=256;
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = [subpaths{mi}, '/' num2str(SystemLength) '/'];  % for reference version
    listr = readmatrix([path, 'list.txt']);
    listrStr = readmatrix([path, 'list.txt'], 'OutputType', 'string');
    listN = listr.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = sprintf('%s',listrStr(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end

    Nsystem = 1;
    %Pair
    % DataIndex = 31; filename='Pair'; Ylab='$N_\mathrm{pairwise} / \mathrm{spin}$'; %% per spin
    % % RND
    % DataIndex = 32; filename='RND'; Ylab='$N_\mathrm{RND}/ \mathrm{spin}$'; %% per spin
    % t_MCS
    DataIndex = 30; filename='t_MCS'; Ylab='$t_{\mathrm{MCS}} / \mathrm{spin}$'; Nsystem = SystemLength ^2; %% per
    % spin
    %Accept
    % DataIndex = 26; filename='Accept'; Ylab='Acceptance Ratio';
    %Accept OR
    % DataIndex = 27; filename='Accept_OR'; Ylab='OR Acceptance Ratio';

    DataToPlot= data(DataIndex,:) / Nsystem; %run/step
    ErrorOfPlot=error(DataIndex,:) / Nsystem;
    % h = errorbar(listr, DataToPlot, ErrorOfPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listr, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 5, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 4
    h.MarkerFaceColor = h.Color;
    % end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    legendStrings = [legendStrings, string([ MethodName])];

end
xlim([0, 35])
% ylim([2e-7 10e-7])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =12;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
ylabel (Ylab, 'FontSize', myFontSize)
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$', 'FontSize', myFontSize)

xlabel('$r_\mathrm{box}$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','northeast','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','north','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','southeast','Interpreter','latex', 'FontSize', 10);
leg1 = legend(legendStrings, 'Location','best','Interpreter','latex', 'FontSize', 10);

set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
% axis 'padded'
% pbaspect([1.5 1 1])
% 
% filename = ['Box_r_Vs_L_clock_', filename];
% filename = ['Box_r_Vs_L_clock_OR_', filename];
% filename = ['Box_r_Vs_L_clock_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_', filename];
% filename = ['Box_r_Vs_L_Diff_fukui_Dythin_', filename];
% filename = ['Box_r_Vs_L_compare_', filename];
filename = ['Box_r_Vs_L_compare_OR_L256_', filename];
exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])

%% set subfigs Fukui-Todo vs Dynamic Thin t_MCS Plain Main
fig = figure;
fig.Position  = [0 0 500 330];
% tiledlayout(4,2)
% axs = tight_subplot(2,2,[.05 .13],[.15 .03],[.13 .03])
axs = tight_subplot(2,2,[.05 .13],[.15 .03],[.12 .03])

clf(axs(1:3))
%%---------------------------------------------------------------------------------
%% Fukui-Todo vs Dynamic Thin t_MCS Plain Main

% clf('reset')
for ii=1:3
% f = figure;
% f.Position = [100 100 300 200];
% methods=["clock", "SCO" , "SCO (OR)" , "Tomita"];
% methods=["32", "64", "128" , "256"];
% SystemLength=[32, 64, 128 , 256];
% methods=["null", "null", "null", "OR-Clock", "null"];

% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];
subpath1 = 'E:\visualize\NewAnalysis\DynamicThinning_Vs_Fukui\Sizes\clock\DyThin';
subpath2 = 'E:\visualize\NewAnalysis\DynamicThinning_Vs_Fukui\Sizes\clock\FT';
subpath3 = 'E:\visualize\NewAnalysis\DynamicThinning_Vs_Fukui\Sizes\SCO\DyThin';
subpath4 = 'E:\visualize\NewAnalysis\DynamicThinning_Vs_Fukui\Sizes\SCO\FT';

% methods=["DyThin", "FT"];
methods = ["Fukui-Todo", "dynamic thinning"];
subpaths={subpath2, subpath1};  %clock
% subpaths={subpath4, subpath3}; %SCO
isOR= false;
cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', '^', 'v'};
legendStrings=[];
tauExponents=[];
mi=0;
SystemLength=256;
axes(axs(ii))
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = [subpaths{mi}, '/' num2str(SystemLength) '/'];  % for reference version
    listr = readmatrix([path, 'list.txt']);
    listrStr = readmatrix([path, 'list.txt'], 'OutputType', 'string');
    listN = listr.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = sprintf('%s',listrStr(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end

    Nsystem = 1;
    if (ii ==1)
        %Pair
        DataIndex = 31; filename='Pair'; Ylab='$N_\mathrm{pw}$'; %% per spin
    elseif (ii == 3)
        % RND
        DataIndex = 32; filename='RND'; Ylab='$N_\mathrm{RND}$'; %% per spin
    elseif (ii == 2)
        % t_MCS
        DataIndex = 30; filename='t_MCS'; Ylab='$t_{\mathrm{spin}} \; [\mu s]$'; Nsystem = SystemLength ^2 * 1e-6; %% per
        % spin
    elseif (ii == 4)
        if (isOR)
            %Accept OR
            DataIndex = 27; filename='Accept_OR'; Ylab='$P_{\mathrm{acc, OR}}$';
        else
            %Accept
            DataIndex = 26; filename='Accept'; Ylab='$P_{\mathrm{acc}}$';
        end
        
    end


    DataToPlot= data(DataIndex,:) / Nsystem; %run/step
    ErrorOfPlot=error(DataIndex,:) / Nsystem;
    % h = errorbar(listr, DataToPlot, ErrorOfPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listr, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 4
    h.MarkerFaceColor = h.Color;
    % end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    legendStrings = [legendStrings, string(MethodName)];

end

% ylim([2e-7 10e-7])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =13;
set(gca,'fontsize',myFontSize) 




% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','northeast','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','north','Interpreter','latex', 'FontSize', 10);
% leg1 = legend(legendStrings, 'Location','southeast','Interpreter','latex', 'FontSize', 10);

if (ii == 1)
    % leg1 = legend(legendStrings, 'Location','southeast','Interpreter','latex', 'FontSize', myFontSize-1);
    % set(leg1,'Box','off')
    set(gca, 'YScale', 'log')   %using plot
    
elseif (ii == 2)
    ylim([0, 1.5])


elseif (ii == 3)
    leg1 = legend(legendStrings, 'Location','best','Interpreter','latex', 'FontSize', myFontSize);
    leg = findobj(gcf, 'Type', 'Legend');
    leg.Position = [0.73 0.16 0.1 0.1];
    set(leg1,'Box','off')



    set(gca, 'YScale', 'log')   %using plot

    yticks([4, 8, 12])
    ylim([4, 16])
     % yticks([1e1, 1e2])
     % ylim([6, 300])
    % if (isOR)
    %     %Accept OR
    % else
    %     yticks([1e1, 1e2])
    % end
elseif (ii == 4)
    if (isOR)
        yticks([ 0.3, 0.6, 0.9])
    else
        %Accept
    end
    
end


% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

ylabel (Ylab, 'FontSize', myFontSize)   % show y label
xlabel('$r_\mathrm{box}$', 'FontSize', myFontSize) % show x label



% set(ax,'xticklabel',[])  % hide x axis  % show y label
% set(ax,'yticklabel',[])  % hide y axis  % show x label


% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot

% pbaspect([1.5 1 1])
% 
end
%%---------------------------------------------------------------------------------
%% remove labels names
set(axs(1),'XTickLabel','');
% % set(axs(2:2:end),'YTickLabel','')
xlabel(axs(1), '');

linkaxes([axs(1:3)],'x')

% add labels to subfigs
ListSubfigLabels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)'}
AddLetters2Plots(fig, ListSubfigLabels, 'HShift', 0.02, 'VShift', 0.0, ...
    'Direction', 'TopDown', 'FontWeight', 'normal', 'FitLocation', true, 'FontSize', myFontSize)



%%---------------------------------------------------------------------------------
%% save

% filename = ['Box_r_Vs_L_clock_OR_', filename];
% filename = ['Box_r_Vs_L_clock_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_', filename];
% filename = ['Box_r_Vs_L_SCO_OR_', filename];

filename = 'Box_r_Vs_L_clock_DyThin_compare_FT_clock';
% filename = 'Box_r_Vs_L_clock_DyThin_compare_FT_SCO';

exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])