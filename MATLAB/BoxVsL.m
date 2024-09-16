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
% fig = figure;
% fig.Position  = [100 100 700 1200];
% tiledlayout(4,2)
ha = tight_subplot(2,2,[.01 .03],[.15 .03],[.15 .03])

%%---------------------------------------------------------------------------------
%% t_MCS Plain Main subfigs

% clf('reset')
ii=4;
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
markerShapes = {'^', 'v', 's', '^', 'v'};
legendStrings=[];
tauExponents=[];
mi=0;
axes(ha(ii))
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
% 

%%---------------------------------------------------------------------------------
%% remove labels names
set(ha(1:end-2),'XTickLabel','');
set(ha(2:2:end),'YTickLabel','')
xlabel(ha(1:end-2), '');
ylabel(ha(2:2:end), '');
% filename = ['Box_r_Vs_L_clock_', filename];
filename = ['Box_r_Vs_L_SCO_', filename];
fig = figure(1);
ListSubfigLabels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)'}
AddLetters2Plots(fig, ListSubfigLabels, 'HShift', 0, 'VShift', 0, ...
    'Direction', 'TopDown', 'FontWeight', 'normal')

% filename = ['Box_r_Vs_L_clock_OR_', filename];
% filename = ['Box_r_Vs_L_clock_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_', filename];
% filename = ['Box_r_Vs_L_SCO_OR_', filename];
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


%%---------------------------------------------------------------------------------
%% Fukui-Todo vs Dynamic Thin t_MCS Plain Main

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

methods=["FT", "DyThin"];
subpaths={subpath3, subpath4};  %clock
% subpaths={subpath10, subpath5}; %SCO

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
    DataIndex = 31; filename='Pair'; Ylab='$C_\mathrm{pairwise}$'; %% per spin
    % % RND
    % DataIndex = 32; filename='RND'; Ylab='$C_\mathrm{RND}$'; %% per spin
    % t_MCS
    % DataIndex = 30; filename='t_MCS'; Ylab='$t_{\mathrm{update}}$'; Nsystem = SystemLength ^2; %% per
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
    legendStrings = [legendStrings, string(MethodName)];

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
set(gca, 'YScale', 'log')   %using plot
% axis 'padded'
% pbaspect([1.5 1 1])
% 
% filename = ['Box_r_Vs_L_clock_', filename];
% filename = ['Box_r_Vs_L_clock_OR_', filename];
% filename = ['Box_r_Vs_L_clock_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_DyThin_', filename];
% filename = ['Box_r_Vs_L_SCO_', filename];
% filename = ['Box_r_Vs_L_Diff_fukui_Dythin_', filename];
filename = ['Box_r_Vs_L_compare_FT_Dythin_clock_', filename];
% filename = ['Box_r_Vs_L_compare_FT_Dythin_SCO_', filename];
exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])