%% Initialize
clear all
clc
addpath(pwd);
addpath([pwd,'\matlab']);
addpath([pwd,'\matlab\mc']);
addpath([pwd,'\matlab\lattice2d']);
addpath([pwd,'\matlab\exportFig']);


date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
save_address = ['E:/visualize/output/', date_time_string , '_Pacc_tMCS'];
[status, msg, msgID] = mkdir(save_address)
dpi = 600;

%% Load locations & auto correlation data
% pathSammple0 = ['../visualize/input/', '0/']

%%---------------------------------------------------------------------------------
%% acceptance

clf('reset')


methods=["Clock", "SCO", "Tomita", "Metropolis"];
legendStrings=[];
for method = methods
    MethodName = sprintf('%s',method);
    path = ['E:/visualize/inputAcceptance/' MethodName '/']
    listL = readmatrix([path, 'list.txt'])
    listN = listL.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = num2str(i);
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end

    plot(listL, data(24,:), 'Marker', '.','MarkerSize', 15);
    hold on
    legendStrings = [legendStrings, method];;
end

ylim([0.2 0.255])
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
set(gca,'fontsize',12) 
ylabel ('Accptance')
xlabel('$L$')
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
hl = legend(legendStrings, 'Location','northeast','Interpreter','latex');

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

exportgraphics(gcf,[save_address '/', 'Acceptance.png'],'Resolution',dpi')

%%---------------------------------------------------------------------------------
%% OR acceptance
clf('reset')


methods=["Clock", "SCO", "Tomita", "Metropolis"];
legendStrings=[];
for method = methods
    MethodName = sprintf('%s',method);
    path = ['E:/visualize/inputAcceptance/' MethodName '/']
    listL = readmatrix([path, 'list.txt'])
    listN = listL.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = num2str(i);
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end

    plot(listL, data(25,:), 'Marker', '.','MarkerSize', 15);
    hold on
    legendStrings = [legendStrings,  string(['OR-' MethodName])];
end

ylim([0.8 1.01])
yl = [0.8 1.00];
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
set(gca,'fontsize',12) 
ylabel ('Accptance')
xlabel('$L$')
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
hl = legend(legendStrings, 'Location','northeast','Interpreter','latex');

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

exportgraphics(gcf,[save_address '/', 'OR Acceptance.png'],'Resolution',dpi)

%%---------------------------------------------------------------------------------
%% t_MCS

clf('reset')


methods=["Metropolis", "Clock", "SCO", "Tomita"];
legendStrings=[];
t_MCSExponents=[];
mi=1;
for method = methods
    MethodName = sprintf('%s',method);
%     path = ['E:/visualize/inputAutoCorrelationTime/' MethodName '/'];
    path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = num2str(i);
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end

%     plot(listL, data(24,:), 'Marker', '.','MarkerSize', 15);
    plot(listL, data(28,:), 'Marker', '.','MarkerSize', 15);
%     errorbar(listL, data(28,:), error(28,:), 'Marker', '.','MarkerSize', 15);
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    xdata = listL(1:N);
    ydata = data(28,1:N);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));
    t_MCSExponents(mi) = p(1);
    mi = mi+1;
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    legendStrings = [legendStrings, string(['$t_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
end

% ylim([0.18 0.27])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
set(gca,'fontsize',12) 
ylabel ('$t_{\mathrm{MCS}}$')
xlabel('$L$')
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
hl = legend(legendStrings, 'Location','southeast','Interpreter','latex');
hl = legend('Location','southeast','Interpreter','latex');

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% set(gcf,'Color',[0 0 0]);
% saveas(gcf, [save_address '/', 'energy.png'])
exportgraphics(gcf,[save_address '/', 't_MCS.png'],'Resolution',dpi)


%%---------------------------------------------------------------------------------
%% t_MCS NON OR

clf('reset')


methods=["Metropolis", "Clock", "SCO", "Tomita"];
legendStrings=[];
t_MCSExponents=[];
mi=1;
for method = methods
    MethodName = sprintf('%s',method);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/']
    path = ['E:/visualize/inputAutoCorrelationTimeNonOR/' MethodName '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = num2str(i);
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end

%     plot(listL, data(24,:), 'Marker', '.','MarkerSize', 15);
    plot(listL, data(28,:), 'Marker', '.','MarkerSize', 15);
%     errorbar(listL, data(28,:), error(28,:), 'Marker', '.','MarkerSize', 15);
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    xdata = listL(1:N);
    ydata = data(28,1:N);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));
    t_MCSExponents(mi) = p(1);
    mi = mi+1;
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    legendStrings = [legendStrings, string(['$t_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
end

% ylim([0.18 0.27])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
set(gca,'fontsize',12) 
ylabel ('$t_{\mathrm{MCS}}$')
xlabel('$L$')
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
hl = legend(legendStrings, 'Location','southeast','Interpreter','latex');
hl = legend('Location','southeast','Interpreter','latex');

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% set(gcf,'Color',[0 0 0]);
% saveas(gcf, [save_address '/', 'energy.png'])
exportgraphics(gcf,[save_address '/', 't_MCS.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'acc.png'],'Resolution',dpi)

%%---------------------------------------------------------------------------------
%% t_MCS  LRON

clf('reset')


methods=["Clock", "SCO", "Tomita"];
legendStrings=[];
for method = methods
    MethodName = sprintf('%s',method);
    path = ['E:/visualize/inputAcceptance/' MethodName '/']
    listL = readmatrix([path, 'list.txt'])
    listN = listL.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = num2str(i);
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end

%     plot(listL, data(24,:), 'Marker', '.','MarkerSize', 15);
    plot(listL, data(28,:)'./(listL.^2.4), 'Marker', '.','MarkerSize', 15);
%     errorbar(listL, data(28,:), error(28,:), 'Marker', '.','MarkerSize', 15);
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
    legendStrings = [legendStrings, method];

end

% ylim([0.18 0.27])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
set(gca,'fontsize',12) 
ylabel ('$t_{\mathrm{MCS}} \; L^{-2.4}$')
xlabel('$L$')
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
hl = legend(legendStrings, 'Location','southeast','Interpreter','latex');
hl = legend('Location','southeast','Interpreter','latex');

% set(hl, 'TextColor','w')
% set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
% set(gcf,'Color',[0 0 0]);
% saveas(gcf, [save_address '/', 'energy.png'])
exportgraphics(gcf,[save_address '/', 'runPerMCStepLRON.png'],'Resolution',dpi)

%%---------------------------------------------------------------------------------
%% acc Main

clf('reset')


% methods=["Clock", "Box-Clock" , "Metropolis", "OR-Clock", "OR-Metropolis"];
% methods=["null", "null", "null", "OR-Clock", "null"];
methods=["clock", "box-clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

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
    path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    % path = ['E:/visualize/inputAutoCorrelationTime/' MethodName '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = num2str(i);
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    
    DataIndex = 24; %Pacc MC
    % DataIndex = 25; %Pacc OR

    DataToPlot=data(DataIndex,:); %run/step

    % h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    if mi >= 4
         h.MarkerFaceColor = h.Color;
    end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 4:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));

    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];;

end

% ylim([10 1e7])
ylim([0 0.26])
yl = ([0 0.25]);
% xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =12;
set(gca,'fontsize',myFontSize) 
% ylabel ('Acceptance rate')
% ylabel ('$t_{\mathrm{MCS}}$')
ylabel ('Acceptance rate', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$')

% xlabel('$L$')
xlabel ('$L$', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
leg1 = legend(legendStrings, 'Location','west','Interpreter','latex', 'FontSize', myFontSize);
set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
% 

exportgraphics(gcf,[save_address '/', 'AccRegular.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/', 'AccRegular.eps'])
% exportgraphics(gcf,[save_address '/', 'AccOR-Clock.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'AccOR-Clock.eps'])


%%---------------------------------------------------------------------------------
%% acc Clock Main

clf('reset')


% methods=["Clock", "Box-Clock" , "Metropolis", "OR-Clock", "OR-Metropolis"];
methods=["null", "null", "null", "OR-clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

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
    path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    % path = ['E:/visualize/inputAutoCorrelationTime/' MethodName '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = num2str(i);
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    
    % DataIndex = 24; %Pacc MC
    DataIndex = 25; %Pacc OR

    DataToPlot=data(DataIndex,:); %run/step

    % h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    if mi >= 4
         h.MarkerFaceColor = h.Color;
    end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 4:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));

    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];;

end

% ylim([10 1e7])
% ylim([0 0.26])
yl = ylim;
% xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =12;
set(gca,'fontsize',myFontSize) 
% ylabel ('Acceptance rate')
% ylabel ('$t_{\mathrm{MCS}}$')
ylabel ('Acceptance rate', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$')

% xlabel('$L$')
xlabel ('$L$', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','west','Interpreter','latex', 'FontSize', myFontSize);
% set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
% 

% exportgraphics(gcf,[save_address '/', 'AccRegular.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'AccRegular.eps'])
exportgraphics(gcf,[save_address '/', 'AccOR-Clock.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/', 'AccOR-Clock.eps'])


%%---------------------------------------------------------------------------------
%% tau Main

clf('reset')


methods=["clock", "box-clock" , "Metropolis", "OR-clock", "OR-Metropolis"];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

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
    path = ['E:/visualize/inputAutoCorrelationTime/' MethodName '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = num2str(i);
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    
    DataIndex = 22; %tau
    % DataIndex = 28; %run/step
    % DataIndex = 24; %Pacc MC
    % DataIndex = 25; %Pacc OR

    DataToPlot=data(DataIndex,:); %run/step

    % h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    if mi >= 4
         h.MarkerFaceColor = h.Color;
    end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 5:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));

    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];;

end

ylim([10 5e6])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =12;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
% ylabel ('$t_{\mathrm{MCS}}$')
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
ylabel ('$\tau$', 'FontSize', myFontSize)

xlabel('$L$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', myFontSize);
set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% 
exportgraphics(gcf,[save_address '/', 'tauOR.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/', 'tauOR.eps'])
% exportgraphics(gcf,[save_address '/', 'runPerMCS.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'runPerMCS.eps'])
% exportgraphics(gcf,[save_address '/', 'AccRegular.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'AccRegular.eps'])
% exportgraphics(gcf,[save_address '/', 'AccOR-Clock.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'AccOR-Clock.eps'])


%%---------------------------------------------------------------------------------
%% t_MCS Main

clf('reset')


methods=["clock", "box-clock" , "Metropolis", "OR-clock", "OR-Metropolis"];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

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
    path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    % path = ['E:/visualize/inputAutoCorrelationTime/' MethodName '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = num2str(i);
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    
    % DataIndex = 22; %tau
    DataIndex = 28; %run/step
    % DataIndex = 24; %Pacc MC
    % DataIndex = 25; %Pacc OR

    DataToPlot=data(DataIndex,:); %run/step

    % h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    if mi >= 4
         h.MarkerFaceColor = h.Color;
    end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 1:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));

    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];;

end

% ylim([10 1e7])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
myFontSize =12;
% set(gca,'fontsize',12) 
set(gca,'fontsize',12) 

% ylabel ('Acceptance rate')
ylabel ('$t_{\mathrm{MCS}} \, (\mathrm{s})$')
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$')

% xlabel('$L$')
xlabel ('$L$', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', myFontSize);
set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% 
% exportgraphics(gcf,[save_address '/', 'tauOR.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'tauOR.eps'])
exportgraphics(gcf,[save_address '/', 'runPerMCS.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/', 'runPerMCS.eps'])
% exportgraphics(gcf,[save_address '/', 'AccRegular.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'AccRegular.eps'])
% exportgraphics(gcf,[save_address '/', 'AccOR-Clock.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'AccOR-Clock.eps'])


%%---------------------------------------------------------------------------------
%% t_eff Main

clf('reset')


methods=["clock", "box-clock" , "Metropolis", "OR-clock", "OR-Metropolis"];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', '^', 'v'};
legendStrings=[];
tauExponents=[];
TauExrapolatedMethod={};
mi=0;
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = ['E:/visualize/inputAutoCorrelationTime/' MethodName '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN);
    data =[];
    error=[];
    for i=1:N
        name = num2str(i);
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    
    DataIndex = 22; %tau
    DataToPlot=data(DataIndex,:);
    N = length(listL);
    range = 4:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));


    % ----------------- here t_eff begins ------------------
    path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    % path = ['E:/visualize/inputAutoCorrelationTime/' MethodName '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = num2str(i);
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    
    % DataIndex = 22; %tau
    DataIndex = 28; %run/step
    % DataIndex = 24; %Pacc MC
    % DataIndex = 25; %Pacc OR

    runPerStepData=data(DataIndex,:); %run/step
    TauExrapolated = exp(polyval(p, log(listL)));
    TauExrapolatedMethod{mi} = TauExrapolated;
    tEffData = 2 .* TauExrapolated .* runPerStepData';
    DataToPlot = tEffData;
    % h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    if mi >= 4
         h.MarkerFaceColor = h.Color;
    end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 4:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    % xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    % [p,S] = polyfit(log(xdata),log(ydata),1);
    % yFit = polyval(p, log(xFit));


    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];;
end

% ylim([10 1e7])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =12;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
% ylabel ('$t_{\mathrm{MCS}}$')
ylabel ('$t_{\mathrm{eff}} \, (\mathrm{s})$', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$')

% xlabel('$L$')
xlabel ('$L$', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', myFontSize);
set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% 
% exportgraphics(gcf,[save_address '/', 'tauOR.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'tauOR.eps'])
% exportgraphics(gcf,[save_address '/', 'runPerMCS.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'runPerMCS.eps'])
% exportgraphics(gcf,[save_address '/', 'AccRegular.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'AccRegular.eps'])
exportgraphics(gcf,[save_address '/', 't_eff.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/', 't_eff.eps'])



%%---------------------------------------------------------------------------------
%% Acceleration Main

clf('reset')


methods=["Clock", "Box-Clock" , "Metropolis", "OR-Clock", "OR-Metropolis"];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', '^', 'v'};
legendStrings=[];
tauExponents=[];
TauExrapolatedMethod={};
teffMethod={};
mi=0;
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);

    % ----------------- tau extrapolation ------------------
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = ['E:/visualize/inputAutoCorrelationTime/' MethodName '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN);
    data =[];
    error=[];
    for i=1:N
        name = num2str(i);
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    
    DataIndex = 22; %tau
    DataToPlot=data(DataIndex,:);
    N = length(listL);
    range = 4:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    ptau = polyfit(log(xdata),log(ydata),1);


    % ----------------- runtime extrapolation ------------------
    path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    % path = ['E:/visualize/inputAutoCorrelationTime/' MethodName '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN);
    data =[];
    error=[];
    for i=1:N
        name = num2str(i);
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    
    DataIndex = 28; %runtime/step
    DataToPlot=data(DataIndex,:);
    N = length(listL);
    range = 4:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    prun = polyfit(log(xdata),log(ydata),1);


    % ----------------- here t_eff begins ------------------

    SizeData=[64, 96, 128, 144, 196, 256, 384, 512, 1024, 2048, 4096, 8192, 2*8192]; %run/step
    listL=SizeData;
    TauExtrpl = exp(polyval(ptau, log(SizeData)));
    runPerStepExtrpl = exp(polyval(prun, log(SizeData)));
    TauExrapolatedMethod{mi} = TauExtrpl;
    tEffData = 2 .* TauExtrpl .* runPerStepExtrpl ./(3600*24);
    % teffMethod{mi}=tEffData;
    teffMethod{mi}=tEffData ./ (3600*24);
    DataToPlot = tEffData ./ (SizeData.^2);
    % h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    % h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    % if strcmp(method, "Over-relaxed Clock")
    %     h.LineWidth = 2.5;
    %     h.MarkerSize = 9;
    % end
    % if mi >= 4
    %      h.MarkerFaceColor = h.Color;
    % end

    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 4:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    % xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    % [p,S] = polyfit(log(xdata),log(ydata),1);
    % yFit = polyval(p, log(xFit));


    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];;
end

acceleration = teffMethod{3} ./ teffMethod{4};
h = plot(listL, acceleration, 'Marker', markerShapes{1} ,'MarkerSize', 7, ...
        'Color',cmap(1,:), 'LineWidth', 0.8);
% ylim([10 1e7])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =12;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
% ylabel ('$t_{\mathrm{MCS}}$')
ylabel ('$t_{\mathrm{eff}}$', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$')

% xlabel('$L$')
xlabel ('$L$', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', myFontSize);
% set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% 
% exportgraphics(gcf,[save_address '/', 'tauOR.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'tauOR.eps'])
% exportgraphics(gcf,[save_address '/', 'runPerMCS.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'runPerMCS.eps'])
% exportgraphics(gcf,[save_address '/', 'AccRegular.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'AccRegular.eps'])
exportgraphics(gcf,[save_address '/', 'acceleration.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/', 'acceleration.eps'])

%---------------------------------------------------------------------------------
%% plot 1-1/a
Rmax = 80;
figure('Units', 'inches', 'Position', [1, 1, 5, 2.5]);
x = linspace(1, Rmax, 1000);  % Generate 1000 equally spaced values from 1 to 100
y = 1 - 1 ./ a;              % Compute the values of the function

plot(x, y, 'LineWidth', 1.5);                  % Plot the function
hold on;

targetX = 20;
[~, pointIndex] = min(abs(x - targetX));
xPoint = x(pointIndex);
yPoint = y(pointIndex);

% Plot a dotted line from the point to the x-axis
plot([xPoint, xPoint], [0, yPoint], '--', 'Color', 'black');

% Plot a dotted line from the point to the y-axis
plot([0, xPoint], [yPoint, yPoint], '--', 'Color', 'black');


xlabel('$r$');
ylabel('$E_r/E_\infty$');
% title('Plot of 1 - 1/a');
% xlim([1 Rmax])

yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/4) :xl(2));
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
set(gca,'fontsize',12) 
% grid off;

% xlabel('a', 'Interpreter', 'latex', 'FontSize', 10);       % Set x-axis label
% ylabel('$1 - \frac{1}{a}$', 'Interpreter', 'latex', 'FontSize', 10);  % Set y-axis label
% title('Plot of $1 - \frac{1}{a}$', 'Interpreter', 'latex', 'FontSize', 10);  % Set title
% grid on;
% box on;
% set(gca, 'FontSize', 10);  % Set font size for axes ticks and labels

exportgraphics(gcf,[save_address '/', '1minus1divr.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/', '1minus1divr.eps'])


% 
% %%---------------------------------------------------------------------------------
% %% different OR Clock Main
% 
% clf('reset')
% 
% 
% methods=["null", "Box-Clock", "1to1", "5to1" , "10to1", "20to1", "50to1", "100to1"];
% MethodNames=["null", "OR = 0", "OR = 1", "OR = 5", "OR = 10", "OR = 20", "OR = 50", "OR = 100"]
% % methods=["null", "null", "null", "OR-Clock", "null"];
% % methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% % methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];
% 
% cmap = linspecer(length(methods));
% markerShapes = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', ...
%            'x', '+', '*', 'o', 'x', 's', 'd', '^', 'v', '>','<', 'p', 'h'};
% % markerShapes = {'0', 'v', 's', '^', 'v', 's','x','o'};
% legendStrings=[];
% tauExponents=[];
% mi=0;
% for method = methods
%     mi = mi+1;
%     if strcmp(method, "null")
%         continue;
%     end
%     % MethodName = sprintf('%s',method);
%     MethodName = MethodNames(mi);
%     % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
%     path = ['E:/visualize/inputAutoCorrelationTimeORClockDifferentRatio/' sprintf('%s',method) '/'];
%     listL = readmatrix([path, 'list.txt']);
%     listN = listL.^2;
%     N = length(listN)
%     data =[];
%     error=[];
%     for i=1:N
%         name = num2str(i);
%         data_path = [path, name , '/'];
% 
%         % [status, msg, msgID] = mkdir(folder_address)
% 
%         data(:,i) = readmatrix([data_path 'MeanX.csv']);
%         error(:,i) = readmatrix([data_path 'errorX.csv']);
%     end
% 
%     DataIndex = 22; %tau
%     % DataIndex = 28; %run/step
%     % DataIndex = 24; %Pacc MC
%     % DataIndex = 25; %Pacc OR
% 
%     DataToPlot=data(DataIndex,:); %run/step
% 
%     % h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
%     %         'Color',cmap(mi,:), 'LineWidth', 0.8);
%     h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
%             'Color',cmap(mi,:), 'LineWidth', 0.8);
%     if strcmp(method, "Over-relaxed Clock")
%         h.LineWidth = 2.5;
%         h.MarkerSize = 9;
%     end
%     if mi >= 5
%          h.MarkerFaceColor = h.Color;
%     end
%     % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
%     hold on
% %     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
% %     legendStrings = [legendStrings, method];
%     N = length(listL);
%     range = 4:N;
%     xdata = listL(range);
%     ydata = DataToPlot(range);
%     xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
%     p = polyfit(log(xdata),log(ydata),1);
%     yFit = polyval(p, log(xFit));
% 
%     mdl = fitlm(log(xdata),log(ydata));
%     tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
%     tauExponents(mi,2) =  mdl.Coefficients.SE(2);
% %     plot(xFit, exp(yFit) , 'LineWidth', 1 );
%     % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
%     % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
%     legendStrings = [legendStrings, string([ MethodName])];;
% 
% end
% 
% x = 24:0.1:60;
% 
% % Calculate the corresponding y values (y = x^2)
% y = 0.12 .* x.^2;
% 
% % Plot the curve with a thick black line
% plot(x, y, 'k--', 'LineWidth', 2);
% % legendStrings = [legendStrings, "$L^{2}$"];;
% 
% xlim([16 64])
% ylim([50 1e4])
% % yl = ylim;
% % xl = xlim;
% % yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% % xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% % set(gca,'fontsize',12) 
% myFontSize =12;
% set(gca,'fontsize',12) 
% 
% % ylabel ('Acceptance rate')
% % ylabel ('$t_{\mathrm{MCS}}$')
% ylabel ('$\tau_M$', 'FontSize', myFontSize)
% % ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% % ylabel ('$\tau$')
% 
% % xlabel('$L$')
% xlabel ('$L$', 'FontSize', 12)
% % legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% % legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', 11);
% set(leg1,'Box','off')
% 
% % set(hl, 'TextColor','w')
% set(gca,'XMinorTick','on','YMinorTick','on');
% ax = gca;
% ax.TickLength(1) = 0.02;
% set(gca,'TickLabelInterpreter','latex');
% set(0,'defaulttextinterpreter','latex');
% 
% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
% % 
% exportgraphics(gcf,[save_address '/', 'tauORRatio.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'tauORRatio.eps'])
% % exportgraphics(gcf,[save_address '/', 'runPerMCS.png'],'Resolution',dpi)
% % exportgraphics(gcf,[save_address '/', 'runPerMCS.eps'])
% % exportgraphics(gcf,[save_address '/', 'AccRegular.png'],'Resolution',dpi)
% % exportgraphics(gcf,[save_address '/', 'AccRegular.eps'])
% % exportgraphics(gcf,[save_address '/', 'AccOR-Clock.png'],'Resolution',dpi)
% % exportgraphics(gcf,[save_address '/', 'AccOR-Clock.eps'])
%% set subfigs different OR Clock Main
fig = figure;
fig.Position  = [0 0 500 330];
% tiledlayout(4,2)
axs = tight_subplot(2,2,[.05 .13],[.15 .03],[.13 .03])


clf(axs(1:3))
% %%---------------------------------------------------------------------------------
%% make plots

methods=["0to1", "1to1", "5to1" , "10to1", "20to1"];
MethodNames=["0", "1", "5", "10", "20",]
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', ...
           'x', '+', '*', 'o', 'x', 's', 'd', '^', 'v', '>','<', 'p', 'h'};
% markerShapes = {'0', 'v', 's', '^', 'v', 's','x','o'};
legendStrings=[];
tauExponents=[];
mi=0;
axes(axs(3))
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    % MethodName = sprintf('%s',method);
    MethodName = MethodNames(mi);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = ['E:\visualize\NewAnalysis\OR Ratio Near group r=11\' sprintf('%s',method) '/'];
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
    
    DataIndex = 24; %tau
    if strcmp(method, "1to0")
        DataIndex = 22; %tau
    end
    % DataIndex = 28; %run/step
    % DataIndex = 24; %Pacc MC
    % DataIndex = 25; %Pacc OR

    DataToPlot=data(DataIndex,:); %run/step

    % h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 5
         h.MarkerFaceColor = h.Color;
    % end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string(['OR / random $=' sprintf('%s',MethodName) , '$' ])];

end

x = 32:0.1:128;
% x = xlim;
% Calculate the corresponding y values (y = x^2)
y = 0.12 .* x.^2;

% Plot the curve with a thick black line
plot(x, y, 'k', 'LineWidth', 2);
% legendStrings = [legendStrings, "$L^{2}$"];;


% ylim([1e2 1e6])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12)
xticks([40,  100]);
% yticks([1e2, 1e4]);
myFontSize =13;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
% ylabel ('$t_{\mathrm{MCS}}$')
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
ylabel ('$\tau \; [\mathrm{MCS}]$', 'FontSize', myFontSize)

xlabel('$L$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', myFontSize);
% set(leg1,'Box','off')
leg1 = legend(legendStrings, 'Location','best','Interpreter','latex', 'FontSize', myFontSize-1);
leg = findobj(gcf, 'Type', 'Legend');
leg.Position = [0.74 0.23 0.1 0.1];
set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% pbaspect([1.2 1 1])
% 
% filename = ['Different_OR_Ratio_clock_', 'tau'];
% exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/',  filename '.eps'])


% 
% %%---------------------------------------------------------------------------------
% %% different OR Clock t_MCS Main

methods=["0to1", "1to1", "5to1" , "10to1", "20to1"];
MethodNames=["0", "1", "5", "10", "20",]
SystemLength=[32, 48, 64, 96, 128];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', ...
           'x', '+', '*', 'o', 'x', 's', 'd', '^', 'v', '>','<', 'p', 'h'};
% markerShapes = {'0', 'v', 's', '^', 'v', 's','x','o'};
legendStrings=[];
tauExponents=[];
mi=0;
axes(axs(1))
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    % MethodName = sprintf('%s',method);
    MethodName = MethodNames(mi);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = ['E:\visualize\NewAnalysis\OR Ratio Near group r=11 t_MCS\' sprintf('%s',method) '/'];
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
    
    Nsystem = 1;
    %Pair
    % DataIndex = 31; filename='Pair'; Ylab='$C_\mathrm{pairwise}$'; %% per spin
    % % RND
    % DataIndex = 32; filename='RND'; Ylab='$C_\mathrm{RND}$'; %% per spin
    % t_MCS
    DataIndex = 30; filename='t_MCS'; Ylab='$t_{\mathrm{MCS}}$'; Nsystem = 1; %% per
    % spin
    %Accept
    % DataIndex = 26; filename='Accept'; Ylab='Acceptance Ratio';
    %Accept OR
    % DataIndex = 27; filename='Accept_OR'; Ylab='OR Acceptance Ratio';

    DataToPlot= data(DataIndex,:) / Nsystem; %run/step
    ErrorOfPlot=error(DataIndex,:) / Nsystem;

    % h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 5
         h.MarkerFaceColor = h.Color;
    % end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string(['OR $=' sprintf('%s',MethodName) , '$' ])];

end


% ylim([1e2 1e6])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
xticks([40,  100]);
yticks([1e-3, 1e-2]);

myFontSize =13;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
ylabel ('$t_{\mathrm{MCS}} \; [\mathrm{s}]$')
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$', 'FontSize', myFontSize)

xlabel('$L$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','southeast','Interpreter','latex', 'FontSize', myFontSize);
% set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% pbaspect([2 1 1])
%
% filename = ['Box_r_Vs_L_Diff_fukui_Dythin_', filename];
% filename = ['Different_OR_Ratio_clock_', filename];
% exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/',  filename '.eps'])


% %%---------------------------------------------------------------------------------
% %% different OR Clock t_eff Main
methods=["0to1", "1to1", "5to1" , "10to1", "20to1"];
MethodNames=["0", "1", "5", "10", "20"]
SystemLength=[32, 48, 64, 96, 128];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', ...
           'x', '+', '*', 'o', 'x', 's', 'd', '^', 'v', '>','<', 'p', 'h'};
% markerShapes = {'0', 'v', 's', '^', 'v', 's','x','o'};
legendStrings=[];
tauExponents=[];
mi=0;
axes(axs(2))
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = MethodNames(mi);
    path = ['E:\visualize\NewAnalysis\OR Ratio Near group r=11\' sprintf('%s',method) '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN)
    data1 =[];
    error1=[];
    for i=1:N
        name = num2str(listL(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data1(:,i) = readmatrix([data_path 'MeanX.csv']);
        error1(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    DataIndex = 24; %tau

    path = ['E:\visualize\NewAnalysis\OR Ratio Near group r=11 t_MCS\' sprintf('%s',method) '/'];
    data2 =[];
    error2=[];
    for i=1:N
        name = num2str(listL(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data2(:,i) = readmatrix([data_path 'MeanX.csv']);
        error2(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    DataIndex = 30; %run/step



    DataToPlot = data1(22,:) .* data2(30,:); %run/step
    % DataToPlot = data1(29,:); %run/step
    % ErrorOfPlot = error1(29,:);
    % h = errorbar(listL, DataToPlot, ErrorOfPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
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
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 1:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));

    mdl = fitlm(log(xdata),log(ydata));
    % xFit=10:10:10000;
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);
    % plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string(['OR $=' sprintf('%s',MethodName) , '$' ])];

end

% ylim([1e2 1e6])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
xticks([40,  100]);
yticks([ 1e0, 1e1, 1e2]);
% set(gca,'fontsize',12) 
myFontSize =13;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
ylabel ('$t_{\mathrm{eff}} \; [\mathrm{s}]$')
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$', 'FontSize', myFontSize)

xlabel('$L$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];


% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% 
% pbaspect([2 1 1])

% newAxSize = [300, 200]; %[width, height] cm
% % Adjust axis and figure sizes
% ax = gca(fig); 
% % ax.Units = 'Centimeters';
% sizeDiff = newAxSize - ax.Position(1:2); 
% % fig.Units = 'centimeters'; 
% fig.Position(1:2) = fig.Position(1:2) + sizeDiff; 
% ax.Position(1:2) = [300, 200];
% movegui(fig)



%% save
set(axs(1),'XTickLabel','');
% % set(axs(2:2:end),'YTickLabel','')
xlabel(axs(1), '');

linkaxes([axs(1:3)],'x')
ListSubfigLabels = {'(a)', '(b)', '(c)'}
AddLetters2Plots(fig, ListSubfigLabels, 'HShift', 0.02, 'VShift', 0, ...
    'Direction', 'TopDown', 'FontWeight', 'normal', 'FontSize', myFontSize)

filename = ['Different_OR_Ratio_clock'];
exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])


%%---------------------------------------------------------------------------------
%% different OR Tomita tau Main

clf('reset')


methods=["0to1", "5to1" , "10to1", "20to1"];
MethodNames=["OR = 0", "OR = 5", "OR = 10", "OR = 20",]
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', ...
           'x', '+', '*', 'o', 'x', 's', 'd', '^', 'v', '>','<', 'p', 'h'};
% markerShapes = {'0', 'v', 's', '^', 'v', 's','x','o'};
legendStrings=[];
tauExponents=[];
mi=0;
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    % MethodName = sprintf('%s',method);
    MethodName = MethodNames(mi);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = ['E:\visualize\NewAnalysis\OR Ratio Near group r=11 Tomita\' sprintf('%s',method) '/'];
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
    
    DataIndex = 24; %tau
    if strcmp(method, "1to0")
        DataIndex = 22; %tau
    end
    % DataIndex = 28; %run/step
    % DataIndex = 24; %Pacc MC
    % DataIndex = 25; %Pacc OR

    DataToPlot=data(DataIndex,:); %run/step

    h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    % h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 5
         h.MarkerFaceColor = h.Color;
    % end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];;

end

x = 32:0.1:128;
% x = xlim;
% Calculate the corresponding y values (y = x^2)
y = 0.12 .* x.^2;

% Plot the curve with a thick black line
plot(x, y, 'k--', 'LineWidth', 2);
% legendStrings = [legendStrings, "$L^{2}$"];;

% xlim([16 64])
% ylim([1e2 1e5])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =12;
set(gca,'fontsize',12) 

% ylabel ('Acceptance rate')
% ylabel ('$t_{\mathrm{MCS}}$')
ylabel ('$\tau$', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$')

% xlabel('$L$')
xlabel ('$L$', 'FontSize', 12)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', 10);
set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
pbaspect([2 1 1])
set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% 
exportgraphics(gcf,[save_address '/', 'tauORRatio.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/', 'tauORRatio.eps'])
% exportgraphics(gcf,[save_address '/', 'runPerMCS.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'runPerMCS.eps'])
% exportgraphics(gcf,[save_address '/', 'AccRegular.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'AccRegular.eps'])
% exportgraphics(gcf,[save_address '/', 'AccOR-Clock.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'AccOR-Clock.eps'])

%---------------------------------------------------------------------------------
%% tau Non OR

clf('reset')


methods=["Metropolis", "Clock", "SCO", "Tomita"];
legendStrings=[];
tauExponents=[];
mi=1;
for method = methods
    MethodName = sprintf('%s',method);
    path = ['E:/visualize/inputAutoCorrelationTimeNonOR/' MethodName '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = num2str(i);
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end

%     plot(listL, data(24,:), 'Marker', '.','MarkerSize', 15);
%     plot(listL, data(22,:), 'Marker', '.','MarkerSize', 15);
    errorbar(listL, data(22,:), error(22,:), 'Marker', '.','MarkerSize', 15);
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 1:N;
    xdata = listL(range);
    ydata = data(22,range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));
    tauExponents(mi) = p(1);
    mi = mi+1;
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
end

x = 0:0.1:listL(end);

% Calculate the corresponding y values (y = x^2)
y = x.^2;

% Plot the curve with a thick black line
% plot(x, y, 'k', 'LineWidth', 2);
xlim([16 64])
% ylim([0.18 0.27])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
set(gca,'fontsize',12) 
ylabel ('$\tau$')
xlabel('$L$')
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
hl = legend(legendStrings, 'Location','southeast','Interpreter','latex');
hl = legend('Location','southeast','Interpreter','latex');

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% set(gcf,'Color',[0 0 0]);
% saveas(gcf, [save_address '/', 'energy.png'])
exportgraphics(gcf,[save_address '/', 'tau Plain.png'],'Resolution',dpi, ...
    'BackgroundColor','none')

%%---------------------------------------------------------------------------------
%% t_eff

clf('reset')


methods=["Metropolis", "Clock", "SCO", "Tomita"];
legendStrings=[];
for method = methods
    MethodName = sprintf('%s',method);
    path = ['E:/visualize/inputAutoCorrelationTime/' MethodName '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = num2str(i);
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    
    myData = data(28,:);
    myError = error(28,:);
    errorbar(listL, myData, myError, 'Marker', '.','MarkerSize', 15);
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 1:N;
    xdata = listL(range);
    ydata = myData(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
end

% ylim([0.18 0.27])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
set(gca,'fontsize',12) 
ylabel ('$\tau$')
xlabel('$L$')
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
hl = legend(legendStrings, 'Location','southeast','Interpreter','latex');
hl = legend('Location','southeast','Interpreter','latex');

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% set(gcf,'Color',[0 0 0]);
% saveas(gcf, [save_address '/', 'energy.png'])
exportgraphics(gcf,[save_address '/', 't_eff.png'],'Resolution',dpi, ...
    'BackgroundColor','none')

%%---------------------------------------------------------------------------------
%% t_eff Exponents by hand

methods=["Metropolis", "Clock", "SCO", "Tomita"];
t_effExponents = t_MCSExponents .* tauExponents;

t_effExponentsTable = array2table(t_effExponents,...
                      'VariableNames',methods)


writetable(t_effExponentsTable,[save_address '/', 't_eff OR Exponents.xlsx']);
% writetable(t_effExponentsTable,[save_address '/', 't_eff Plain Exponents.xlsx']);


%%---------------------------------------------------------------------------------
%% t_eff simple extrapolate

clf('reset')


methods=["clock", "box-clock" , "Metropolis", "OR-clock", "OR-Metropolis"];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', '^', 'v'};
legendStrings=[];
tauExponents=[];
TauExrapolatedMethod={};
mi=0;
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = ['E:/visualize/inputAutoCorrelationTime/' MethodName '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN);
    data =[];
    error=[];
    for i=1:N
        name = num2str(i);
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    
    DataIndex = 22; %tau
    DataToPlot=data(DataIndex,:);
    N = length(listL);
    range = 4:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));


    % ----------------- here t_eff begins ------------------
    path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    % path = ['E:/visualize/inputAutoCorrelationTime/' MethodName '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN)
    data =[];
    error=[];
    for i=1:N
        name = num2str(i);
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data(:,i) = readmatrix([data_path 'MeanX.csv']);
        error(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    
    % DataIndex = 22; %tau
    DataIndex = 28; %run/step
    % DataIndex = 24; %Pacc MC
    % DataIndex = 25; %Pacc OR

    runPerStepData=data(DataIndex,:); %run/step
    TauExrapolated = exp(polyval(p, log(listL)));
    TauExrapolatedMethod{mi} = TauExrapolated;
    tEffData = 2 .* TauExrapolated .* runPerStepData';
    DataToPlot = tEffData;
    % h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    % h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    % if strcmp(method, "Over-relaxed Clock")
    %     h.LineWidth = 2.5;
    %     h.MarkerSize = 9;
    % end
    % if mi >= 4
    %      h.MarkerFaceColor = h.Color;
    % end
    % hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 4:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    [p,S] = polyfit(log(xdata),log(ydata),1);

    xFit = 10:10:1000;
    yFit = polyval(p, log(xFit));


    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);
    plot(xFit, exp(yFit) , 'LineWidth', 1 );
    hold on;
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];
end

% ylim([10 1e7])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
myFontSize =12;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
% ylabel ('$t_{\mathrm{MCS}}$')
ylabel ('$t_{\mathrm{eff}} \, (\mathrm{s})$', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$')

% xlabel('$L$')
xlabel ('$L$', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', myFontSize);
set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot



% 
% exportgraphics(gcf,[save_address '/', 'tauOR.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'tauOR.eps'])
% exportgraphics(gcf,[save_address '/', 'runPerMCS.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'runPerMCS.eps'])
% exportgraphics(gcf,[save_address '/', 'AccRegular.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'AccRegular.eps'])
exportgraphics(gcf,[save_address '/', 't_eff.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/', 't_eff.eps'])


%% set subfigs compare methods initial Main
fig = figure;
fig.Position  = [0 0 500 330];
% tiledlayout(4,2)
axs = tight_subplot(2,2,[.05 .15],[.15 .03],[.14 .02])


clf(axs(1:3))
%%---------------------------------------------------------------------------------
%% tau initial Main
% Set the paper size to A4[210mm x 297mm]
% a4_width = 21.0; % inches
% a4_height = 29.7; % inches
% figure('Units', 'centimeters', 'Position', [0, 0, 9, 4]);
% set(gcf, 'PaperSize', [a4_width, a4_height]);
% set(gcf, 'PaperPositionMode', 'auto');
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPosition', [0, 0, a4_width, a4_height]);

% clf('reset')

% fig = figure;
% fig.Position = [100 100 400 300];

% methods=["clock", "SCO" , "SCO (OR)", "Tomita", "Metropolis"];
methods=["clock", "SCO" , "Tomita", "Metropolis", "Metropolis (OR)"]
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', 'o', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};
legendStrings=[];
tauExponents=[];
mi=0;
axes(axs(3))
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = ['E:/visualize/NewAnalysis/Tau Different Algorithms/' MethodName '/'];
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
    
    if strcmp(method, "Metropolis (OR)")
        DataIndex = 24; %tau
    else
        DataIndex = 22; %tauM
    end
    
    % DataIndex = 28; %run/step
    % DataIndex = 24; %Pacc MC
    % DataIndex = 25; %Pacc OR

    DataToPlot=data(DataIndex,:); %run/step

    h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    % h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
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
    N = length(listL);
    range = 3:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));

    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];;

end

% ylim([1e2 1e6])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
xticks([20, 40, 90]);
yticks([1e2, 1e4, 1e6]);

myFontSize =13;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
% ylabel ('$t_{\mathrm{MCS}}$')
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
ylabel ('$\tau \; [\mathrm{MCS}]$', 'FontSize', myFontSize)

xlabel('$L$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', myFontSize);
% set(leg1,'Box','off')
leg1 = legend(legendStrings, 'Location','best','Interpreter','latex', 'FontSize', myFontSize-1);
leg = findobj(gcf, 'Type', 'Legend');
leg.Position = [0.74 0.23 0.1 0.1];
set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% pbaspect([1.2 1 1])
% 
% exportgraphics(gcf,[save_address '/', 'PlainTau.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'PlainTau.eps'])

%%---------------------------------------------------------------------------------
%% t_MCS initial Main

% clf('reset')
% fig = figure;
% fig.Position = [100 100 400 300]
% methods=["clock", "SCO" , "SCO (OR)", "Tomita", "Metropolis"];
methods=["clock", "SCO" , "Tomita", "Metropolis", "Metropolis (OR)"]
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', 'o', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};
legendStrings=[];
tauExponents=[];
mi=0;
axes(axs(1))
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = ['E:\visualize\NewAnalysis\t_MCS\different sizes/' MethodName '/'];
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
    
    DataIndex = 30; %run/step
    % DataIndex = 24; %Pacc MC
    % DataIndex = 25; %Pacc OR

    DataToPlot=data(DataIndex,:); %run/step

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
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 3:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));

    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];;

end

% ylim([1e2 1e6])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
xticks([20, 40, 90]);
yticks([1e-4, 1e-2]);

myFontSize =13;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
ylabel ('$t_{\mathrm{MCS}} \; [\mathrm{s}]$')
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$', 'FontSize', myFontSize)

xlabel('$L$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','southeast','Interpreter','latex', 'FontSize', myFontSize);
% set(leg1,'Box','off')


% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% pbaspect([2 1 1])


%%---------------------------------------------------------------------------------
%% t_eff initial Main

% clf('reset')
% fig = figure;
% fig.Position = [100 100 640 300]
% methods=["clock", "SCO" , "SCO (OR)", "Tomita", "Metropolis"];
methods=["clock", "SCO" , "Tomita", "Metropolis", "Metropolis (OR)"]

% methods=["clock", "SCO" , "Tomita"];
% methods=[ "SCO" , "Tomita"];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', 'o', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};legendStrings=[];
tauExponents=[];
mi=0;
axes(axs(2))
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    path = ['E:/visualize/NewAnalysis/Tau Different Algorithms/' MethodName '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN)
    data1 =[];
    error1=[];
    for i=1:N
        name = num2str(listL(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data1(:,i) = readmatrix([data_path 'MeanX.csv']);
        error1(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    DataIndex = 22; %tau

    path = ['E:\visualize\NewAnalysis\t_MCS\different sizes/' MethodName '/'];
    data2 =[];
    error2=[];
    for i=1:N
        name = num2str(listL(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data2(:,i) = readmatrix([data_path 'MeanX.csv']);
        error2(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    DataIndex = 30; %run/step
    

    
    if strcmp(method, "Metropolis (OR)")
        DataIndex1 = 24; %tau
    else
        DataIndex1 = 22; %tauM
    end

    DataToPlot = 2 * data1(DataIndex1,:) .* data2(30,:); %run/step
    % DataToPlot = data1(29,:); %run/step
    % ErrorOfPlot = error1(29,:);
    % h = errorbar(listL, DataToPlot, ErrorOfPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
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
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 3:N;
    if strcmp(method, "Metropolis (OR)")
        range = 3:N;
    end
    xdata = listL(range);
    ydata = DataToPlot(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));

    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);

    xFit = 10:10:5000;
    yFit = polyval(p, log(xFit));
    % plot(xFit, exp(yFit) , 'LineWidth', 1 );
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];;

end

ylim([1e-2 1e3])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
xticks([20, 40, 90]);
yticks([ 1e-2, 1e0, 1e2]);
% set(gca,'fontsize',12) 
myFontSize =13;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
ylabel ('$t_{\mathrm{eff}} \; [\mathrm{s}]$')
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$', 'FontSize', myFontSize)

xlabel('$L$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','best','Interpreter','latex', 'FontSize', myFontSize-1);
% leg = findobj(gcf, 'Type', 'Legend');
% leg.Position = [0.62 0.25 0.1 0.1];
% set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot




%% save
set(axs(1),'XTickLabel','');
% % set(axs(2:2:end),'YTickLabel','')
xlabel(axs(1), '');

linkaxes([axs(1:3)],'x')
ListSubfigLabels = {'(a)', '(b)', '(c)'}
AddLetters2Plots(fig, ListSubfigLabels, 'HShift', 0.02, 'VShift', 0, ...
    'Direction', 'TopDown', 'FontWeight', 'normal', 'FontSize', myFontSize)

filename = ['Compare_Methods_Plain'];
exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])



%% set subfigs compare methods Final Main
fig = figure;
fig.Position  = [0 0 1000 200];
% tiledlayout(4,2)
axs = tight_subplot(1,3,[.1 .07],[.24 .03],[.07 .2])


% clf(axs(1:3))
%%---------------------------------------------------------------------------------
%% different methods final tau Main

% fig = figure;
% fig.Position = [100 100 400 300]
% methods=["0to1", "1to1", "5to1" , "10to1", "20to1"];
% MethodNames=["0", "1", "5", "10", "20",]
% methods=["null", "null", "null", "OR-Clock", "null"];
methods=["clock (box+OR)", "SCO (box+OR)" , "Tomita (box+OR)", "null" , "Metropolis (OR)"];
MethodNames=methods;
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', 'o', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};
legendStrings=[];
tauExponents=[];
mi=0;
axes(axs(2))
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    % MethodName = sprintf('%s',method);
    MethodName = MethodNames(mi);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = ['E:\visualize\NewAnalysis\Different Methods OR final tau\' sprintf('%s',method) '/'];
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
    
    DataIndex = 24; %tau
    % if strcmp(method, "1to0")
    %     DataIndex = 22; %tau
    % end
    % DataIndex = 28; %run/step
    % DataIndex = 24; %Pacc MC
    % DataIndex = 25; %Pacc OR

    DataToPlot=data(DataIndex,:); %run/step

    % h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);

     h.MarkerFaceColor = h.Color;
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on

    range = 1:N;
    xdata = listL(range);
    ydata = DataToPlot(range);

    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);

%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];;

end


% ylim([1e2 1e6])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12)
xticks([40,  100]);
yticks([1e2, 1e3]);
myFontSize =13;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
% ylabel ('$t_{\mathrm{MCS}}$')
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
ylabel ('$\tau \; [\mathrm{MCS}]$', 'FontSize', myFontSize)

xlabel('$L$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', myFontSize);
% set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% pbaspect([1.2 1 1])
% 
% ListSubfigLabels = {'(b)'};
% AddLetters2Plots(fig, ListSubfigLabels, 'HShift', 0.02, 'VShift', 0, ...
%     'Direction', 'TopDown', 'FontWeight', 'normal', 'FontSize', myFontSize)

% filename = ['Different_Method_final_', 'tau'];
% exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/',  filename '.eps'])




%%---------------------------------------------------------------------------------
%% different methods final t_MCS Main
% clf('reset')

% fig = figure;
% fig.Position = [100 100 400 300]
methods=["clock (box+OR)", "SCO (box+OR)" , "Tomita (box+OR)", "null" , "Metropolis (OR)"];
MethodNames=methods;
% SystemLength=[32, 48, 64, 96, 128];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', 'o', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};
legendStrings=[];
tauExponents=[];
mi=0;
Li = 1;
axes(axs(1))

for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    % MethodName = sprintf('%s',method);
    MethodName = MethodNames(mi);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = ['E:\visualize\NewAnalysis\Different Methods OR final t_MCS\' sprintf('%s',method) '/'];
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
    
    Nsystem = 1;
    %Pair
    % DataIndex = 31; filename='Pair'; Ylab='$N_\mathrm{pairwise} / \mathrm{spin}$'; %% per spin
    % % RND
    % DataIndex = 32; filename='RND'; Ylab='$N_\mathrm{RND}/ \mathrm{spin}$'; %% per spin
    % t_MCS
    % DataIndex = 30; filename='t_MCS'; Ylab='$t_{\mathrm{MCS}} / \mathrm{spin}$'; Nsystem = SystemLength(mi) ^2; %% per
    % spin
     DataIndex = 30; filename='t_MCS'; Ylab='$t_{\mathrm{MCS}}$';
    %Accept
    % DataIndex = 26; filename='Accept'; Ylab='Acceptance Ratio';
    %Accept OR
    % DataIndex = 27; filename='Accept_OR'; Ylab='OR Acceptance Ratio';

    DataToPlot= data(DataIndex,:) / Nsystem; %run/step
    ErrorOfPlot=error(DataIndex,:) / Nsystem;
    % 
    % h = errorbar(listL, data(DataIndex,:), error(DataIndex,:), 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 5
         h.MarkerFaceColor = h.Color;
    % end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on

    range = 1:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    
    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);

%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];;

end

% ylim([1e2 1e6])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
xticks([40,  100]);
yticks([1e-3, 1e-2]);

myFontSize =13;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
ylabel ('$t_{\mathrm{MCS}} \; [\mathrm{s}]$')
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$', 'FontSize', myFontSize)

xlabel('$L$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','southeast','Interpreter','latex', 'FontSize', myFontSize);
% set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot

% pbaspect([2 1 1])
%
% ListSubfigLabels = {'(a)'}
% AddLetters2Plots(fig, ListSubfigLabels, 'HShift', 0.02, 'VShift', 0, ...
%     'Direction', 'TopDown', 'FontWeight', 'normal', 'FontSize', myFontSize)


% filename = ['Different_Method_final_', filename];
% exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/',  filename '.eps'])



%%---------------------------------------------------------------------------------
%% different methods final t_eff Main

% fig = figure;
% fig.Position = [100 100 650 300]
methods=["clock (box+OR)", "SCO (box+OR)" , "Tomita (box+OR)", "null" , "Metropolis (OR)"];
MethodNames=methods;
SystemLength=[32, 48, 64, 96, 128];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', 'o', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};
legendStrings=[];
tauExponents=[];
mi=0;
axes(axs(3))
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = MethodNames(mi);
    path = ['E:\visualize\NewAnalysis\Different Methods OR final tau\' sprintf('%s',method) '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN)
    data1 =[];
    error1=[];
    for i=1:N
        name = num2str(listL(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data1(:,i) = readmatrix([data_path 'MeanX.csv']);
        error1(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    DataIndex = 24; %tau

    path = ['E:\visualize\NewAnalysis\Different Methods OR final t_MCS\' sprintf('%s',method) '/'];
    data2 =[];
    error2=[];
    for i=1:N
        name = num2str(listL(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data2(:,i) = readmatrix([data_path 'MeanX.csv']);
        error2(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    DataIndex = 30; %run/step



    DataToPlot = 2 * data1(22,:) .* data2(30,:); %run/step
    % DataToPlot = data1(29,:); %run/step
    % ErrorOfPlot = error1(29,:);
    % h = errorbar(listL, DataToPlot, ErrorOfPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
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
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 1:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));


    xFit = 10:10:5000;
    yFit = polyval(p, log(xFit));
    % plot(xFit, exp(yFit) , 'LineWidth', 1 );

    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];

end

ylim([0.2 200])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
xticks([40,  100]);

yticks([ 1e0, 1e2]);
% set(gca,'fontsize',12) 
myFontSize =13;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
ylabel ('$t_{\mathrm{eff}} \; [\mathrm{s}]$')
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$', 'FontSize', myFontSize)

xlabel('$L$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
leg1 = legend(legendStrings, 'Location','best','Interpreter','latex', 'FontSize', myFontSize);
leg = findobj(gcf, 'Type', 'Legend');
leg.Position = [0.84 0.39 0.1 0.1];
set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% 
% ListSubfigLabels = {'(c)'};
% AddLetters2Plots(fig, ListSubfigLabels, 'HShift', 0.02, 'VShift', 0, ...
%     'Direction', 'TopDown', 'FontWeight', 'normal', 'FontSize', myFontSize)

% filename = ['Different_Method_final_', 't_eff'];
% exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/',  filename '.eps'])


%% save
% set(axs(1),'XTickLabel','');
% % set(axs(2:2:end),'YTickLabel','')
% xlabel(axs(1), '');

linkaxes([axs(1:3)],'x')
ListSubfigLabels = {'(a)', '(b)', '(c)'}
AddLetters2Plots(fig, ListSubfigLabels, 'HShift', 0.005, 'VShift', 0, ...
    'Direction', 'LeftRight', 'FontWeight', 'normal', 'FontSize', myFontSize)

filename = ['Different_Method_final'];
exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])

%% set subfigs compare methods Final 2
fig = figure;
fig.Position  = [0 0 500 330];
% tiledlayout(4,2)
axs = tight_subplot(2,2,[.05 .15],[.15 .03],[.14 .02])


clf(axs(1:3))
%%---------------------------------------------------------------------------------
%% tau Final 2
% Set the paper size to A4[210mm x 297mm]
% a4_width = 21.0; % inches
% a4_height = 29.7; % inches
% figure('Units', 'centimeters', 'Position', [0, 0, 9, 4]);
% set(gcf, 'PaperSize', [a4_width, a4_height]);
% set(gcf, 'PaperPositionMode', 'auto');
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPosition', [0, 0, a4_width, a4_height]);

% clf('reset')

% fig = figure;
% fig.Position = [100 100 400 300];

% methods=["clock", "SCO" , "SCO (OR)", "Tomita", "Metropolis"];
methods=["clock (box+OR)", "SCO (box+OR)" , "Tomita (box+OR)", "null" , "Metropolis (OR)"];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', 'o', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};
legendStrings=[];
tauExponents=[];
mi=0;
axes(axs(3))
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = ['E:\visualize\NewAnalysis\Different Methods OR final tau\' sprintf('%s',method) '/'];
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
    
    if strcmp(method, "Metropolis (OR)")
        DataIndex = 24; %tau
    else
        DataIndex = 22; %tauM
    end
    
    % DataIndex = 28; %run/step
    % DataIndex = 24; %Pacc MC
    % DataIndex = 25; %Pacc OR

    DataToPlot=data(DataIndex,:); %run/step

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
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 1:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));

    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];;

end

% ylim([1e2 1e6])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
xticks([40,  100]);
yticks([1e2, 1e3]);

myFontSize =13;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
% ylabel ('$t_{\mathrm{MCS}}$')
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
ylabel ('$\tau \; [\mathrm{MCS}]$', 'FontSize', myFontSize)

xlabel('$L$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','northwest','Interpreter','latex', 'FontSize', myFontSize);
% set(leg1,'Box','off')
leg1 = legend(legendStrings, 'Location','best','Interpreter','latex', 'FontSize', myFontSize);
leg = findobj(gcf, 'Type', 'Legend');
leg.Position = [0.74 0.23 0.1 0.1];
set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% pbaspect([1.2 1 1])
% 
% exportgraphics(gcf,[save_address '/', 'PlainTau.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'PlainTau.eps'])

%%---------------------------------------------------------------------------------
%% t_MCS Final 2

% clf('reset')
% fig = figure;
% fig.Position = [100 100 400 300]
% methods=["clock", "SCO" , "SCO (OR)", "Tomita", "Metropolis"];
methods=["clock (box+OR)", "SCO (box+OR)" , "Tomita (box+OR)", "null" , "Metropolis (OR)"];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', 'o', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};
legendStrings=[];
tauExponents=[];
mi=0;
axes(axs(1))
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    % path = ['E:/visualize/inputAcceptance/' MethodName '/'];
    path = ['E:\visualize\NewAnalysis\Different Methods OR final t_MCS\' sprintf('%s',method) '/'];
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
    
    DataIndex = 30; %run/step
    % DataIndex = 24; %Pacc MC
    % DataIndex = 25; %Pacc OR

    DataToPlot=data(DataIndex,:); %run/step

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
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 1:N;
    xdata = listL(range);
    ydata = DataToPlot(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));

    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];;

end

% ylim([1e2 1e6])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% set(gca,'fontsize',12) 
xticks([40,  100]);
yticks([1e-3, 1e-2]);

myFontSize =13;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
ylabel ('$t_{\mathrm{MCS}} \; [\mathrm{s}]$')
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$', 'FontSize', myFontSize)

xlabel('$L$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','southeast','Interpreter','latex', 'FontSize', myFontSize);
% set(leg1,'Box','off')


% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
% pbaspect([2 1 1])


%%---------------------------------------------------------------------------------
%% t_eff Final 2

% clf('reset')
% fig = figure;
% fig.Position = [100 100 640 300]
% methods=["clock", "SCO" , "SCO (OR)", "Tomita", "Metropolis"];
methods=["clock (box+OR)", "SCO (box+OR)" , "Tomita (box+OR)", "null" , "Metropolis (OR)"];

% methods=["clock", "SCO" , "Tomita"];
% methods=[ "SCO" , "Tomita"];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', 'o', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};legendStrings=[];
tauExponents=[];
mi=0;
axes(axs(2))
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    path = ['E:\visualize\NewAnalysis\Different Methods OR final tau\' sprintf('%s',method) '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN)
    data1 =[];
    error1=[];
    for i=1:N
        name = num2str(listL(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data1(:,i) = readmatrix([data_path 'MeanX.csv']);
        error1(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    DataIndex = 22; %tau

    path = ['E:\visualize\NewAnalysis\Different Methods OR final t_MCS\' sprintf('%s',method) '/'];
    data2 =[];
    error2=[];
    for i=1:N
        name = num2str(listL(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data2(:,i) = readmatrix([data_path 'MeanX.csv']);
        error2(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    DataIndex = 30; %run/step
    

    
    if strcmp(method, "Metropolis (OR)")
        DataIndex1 = 24; %tau
    else
        DataIndex1 = 22; %tauM
    end

    DataToPlot = 2 * data1(DataIndex1,:) .* data2(30,:); %run/step
    % DataToPlot = data1(29,:); %run/step
    % ErrorOfPlot = error1(29,:);
    % h = errorbar(listL, DataToPlot, ErrorOfPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
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
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 1:N;
    if strcmp(method, "Metropolis (OR)")
        range = 3:N;
    end
    xdata = listL(range);
    ydata = DataToPlot(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));

    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);

    xFit = 10:10:5000;
    yFit = polyval(p, log(xFit));
    % plot(xFit, exp(yFit) , 'LineWidth', 1 );
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];;

end

ylim([1e-2 1e3])
% ylim([0 0.3])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
xticks([40,  100]);

yticks([ 1e0, 1e2]);
% set(gca,'fontsize',12) 
myFontSize =13;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
ylabel ('$t_{\mathrm{eff}} \; [\mathrm{s}]$')
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$', 'FontSize', myFontSize)

xlabel('$L$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% leg1 = legend(legendStrings, 'Location','best','Interpreter','latex', 'FontSize', myFontSize-1);
% leg = findobj(gcf, 'Type', 'Legend');
% leg.Position = [0.62 0.25 0.1 0.1];
% set(leg1,'Box','off')

% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot




%% save
set(axs(1),'XTickLabel','');
% % set(axs(2:2:end),'YTickLabel','')
xlabel(axs(1), '');

linkaxes([axs(1:3)],'x')
ListSubfigLabels = {'(a)', '(b)', '(c)'}
AddLetters2Plots(fig, ListSubfigLabels, 'HShift', 0.02, 'VShift', 0, ...
    'Direction', 'TopDown', 'FontWeight', 'normal', 'FontSize', myFontSize)

filename = ['Different_Method_final'];
exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])

%% t_eff plain vs final
fig = figure;
fig.Position  = [0 0 500 390];
% tiledlayout(4,2)
% axs = tight_subplot(2,2,[.05 .15],[.15 .03],[.14 .02])

% % clf(axs(1:3))
% %%---------------------------------------------------------------------------------
% %% t_eff Plain in plain vs final

% clf('reset')
% fig = figure;
% fig.Position = [100 100 640 300]
% methods=["clock", "SCO" , "SCO (OR)", "Tomita", "Metropolis"];
methods=["clock", "SCO" , "Tomita"]

% methods=["clock", "SCO" , "Tomita"];
% methods=[ "SCO" , "Tomita"];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'--^', '--v', '--s', 'o', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};legendStrings=[];
tauExponents=[];
mi=0;
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    path = ['E:/visualize/NewAnalysis/Tau Different Algorithms/' MethodName '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN)
    data1 =[];
    error1=[];
    for i=1:N
        name = num2str(listL(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data1(:,i) = readmatrix([data_path 'MeanX.csv']);
        error1(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    DataIndex = 22; %tau

    path = ['E:\visualize\NewAnalysis\t_MCS\different sizes/' MethodName '/'];
    data2 =[];
    error2=[];
    for i=1:N
        name = num2str(listL(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data2(:,i) = readmatrix([data_path 'MeanX.csv']);
        error2(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    DataIndex = 30; %run/step
    

    
    if strcmp(method, "Metropolis (OR)")
        DataIndex1 = 24; %tau
    else
        DataIndex1 = 22; %tauM
    end

    DataToPlot = 2 * data1(DataIndex1,:) .* data2(30,:); %run/step
    % DataToPlot = data1(29,:); %run/step
    % ErrorOfPlot = error1(29,:);
    % h = errorbar(listL, DataToPlot, ErrorOfPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listL, DataToPlot,  markerShapes{mi} ,'MarkerSize', 7, ...
            'Color',cmap(mi,:), 'LineWidth', 0.8);
    if strcmp(method, "Over-relaxed Clock")
        h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 4
    % h.MarkerFaceColor = h.Color;
    % end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 1:N;
    if strcmp(method, "Metropolis (OR)")
        range = 3:N;
    end
    xdata = listL(range);
    ydata = DataToPlot(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));

    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);

    xFit = 10:10:5000;
    yFit = polyval(p, log(xFit));
    % plot(xFit, exp(yFit) , 'LineWidth', 1 );
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];;

end
% %%---------------------------------------------------------------------------------
% %% t_eff Final in plain vs final

% clf('reset')
% fig = figure;
% fig.Position = [100 100 640 300]
% methods=["clock", "SCO" , "SCO (OR)", "Tomita", "Metropolis"];
methods=["clock (box+OR)", "SCO (box+OR)" , "Tomita (box+OR)", "null"];

% methods=["clock", "SCO" , "Tomita"];
% methods=[ "SCO" , "Tomita"];
% methods=["null", "null", "null", "OR-Clock", "null"];
% methods=["Clock", "Box-Clock" , "Metropolis","null", "null"];
% methods=[ "Metropolis", "Clock", "Over-relaxation", "null" , "Over-relaxed Clock"];

cmap = linspecer(length(methods));
markerShapes = {'^', 'v', 's', 'o', 'd', '>', '<', 'p', 'h', 'x', '+', '*'};
tauExponents=[];
mi=0;
for method = methods
    mi = mi+1;
    if strcmp(method, "null")
        continue;
    end
    MethodName = sprintf('%s',method);
    path = ['E:\visualize\NewAnalysis\Different Methods OR final tau\' sprintf('%s',method) '/'];
    listL = readmatrix([path, 'list.txt']);
    listN = listL.^2;
    N = length(listN)
    data1 =[];
    error1=[];
    for i=1:N
        name = num2str(listL(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data1(:,i) = readmatrix([data_path 'MeanX.csv']);
        error1(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    DataIndex = 22; %tau

    path = ['E:\visualize\NewAnalysis\Different Methods OR final t_MCS\' sprintf('%s',method) '/'];
    data2 =[];
    error2=[];
    for i=1:N
        name = num2str(listL(i));
        data_path = [path, name , '/'];

        % [status, msg, msgID] = mkdir(folder_address)

        data2(:,i) = readmatrix([data_path 'MeanX.csv']);
        error2(:,i) = readmatrix([data_path 'errorX.csv']);
    end
    DataIndex = 30; %run/step
    

    
    if strcmp(method, "Metropolis (OR)")
        DataIndex1 = 24; %tau
    else
        DataIndex1 = 22; %tauM
    end

    DataToPlot = 2 * data1(DataIndex1,:) .* data2(30,:); %run/step
    % DataToPlot = data1(29,:); %run/step
    % ErrorOfPlot = error1(29,:);
    % h = errorbar(listL, DataToPlot, ErrorOfPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 7, ...
    %         'Color',cmap(mi,:), 'LineWidth', 0.8);
    h = plot(listL, DataToPlot, 'Marker', markerShapes{mi} ,'MarkerSize', 8, ...
            'Color',cmap(mi,:), 'LineWidth', 1.1);

    if strcmp(method, "clock (box+OR)")
        % h.LineWidth = 2.5;
        h.MarkerSize = 9;
    end
    % if mi >= 4
    h.MarkerFaceColor = h.Color;
    % end
    % set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    hold on
%     legendStrings = [legendStrings, string(['$P_{\mbox{' MethodName '}}$'])];;
%     legendStrings = [legendStrings, method];
    N = length(listL);
    range = 1:N;
    if strcmp(method, "Metropolis (OR)")
        range = 3:N;
    end
    xdata = listL(range);
    ydata = DataToPlot(range);
    xFit = linspace(min(xdata) - min(xdata)/10, max(xdata) + max(xdata)/10, 500);
    p = polyfit(log(xdata),log(ydata),1);
    yFit = polyval(p, log(xFit));

    mdl = fitlm(log(xdata),log(ydata));
    tauExponents(mi,1) = mdl.Coefficients.Estimate(2);
    tauExponents(mi,2) =  mdl.Coefficients.SE(2);

    xFit = 10:10:5000;
    yFit = polyval(p, log(xFit));
    % plot(xFit, exp(yFit) , 'LineWidth', 1 );
%     plot(xFit, exp(yFit) , 'LineWidth', 1 );
    % legendStrings = [legendStrings, string(['$\tau_{\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];
    % legendStrings = [legendStrings, string(['${\mathrm{' MethodName '}} \sim L^{ ' num2str(p(1), 3) '}$'])];;
    legendStrings = [legendStrings, string([ MethodName])];;

end
legendStrings = ["clock", "SCO" , "Tomita" , ...
                 "clock (box+OR)", "SCO (box+OR)" , "Tomita (box+OR)"];
xlim([10 190])
ylim([7e-2 2e4])
% yl = ylim;
% xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
% xticks([]);

yticks([1e-1, 1e0, 1e1, 1e2, 1e3]);
% xticks([20 , 100]);
% set(gca,'fontsize',12) 
myFontSize =13;
set(gca,'fontsize',myFontSize) 

% ylabel ('Acceptance rate')
ylabel ('$t_{\mathrm{eff}} \; [\mathrm{s}]$')
% ylabel ('Microstate generation time', 'FontSize', myFontSize)
% ylabel ('Correlation between microstates', 'FontSize', myFontSize)
% ylabel ('$\tau$', 'FontSize', myFontSize)

xlabel('$L$', 'FontSize', myFontSize)
% xlabel ('Size of System', 'FontSize', myFontSize)
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
leg1 = legend(legendStrings, 'Location','northwest','Orientation', 'vertical', ...
             'NumColumns', 2, 'Interpreter','latex', 'FontSize', myFontSize-2);
% leg = findobj(gcf, 'Type', 'Legend');
% leg.Position = [0.62 0.25 0.1 0.1];
set(leg1,'Box','off')
lgd1.NumColumns=2;

% set(hl, 'TextColor','w')
set(gca, 'XScale', 'log')   %using plot
set(gca, 'YScale', 'log')   %using plot
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.03;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
%% save
filename = ['t_eff_Final_Vs_Plain'];
exportgraphics(gcf,[save_address '/',  filename '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  filename '.eps'])


%%---------------------------------------------------------------------------------
%% constant box Vs near neighbors

clf('reset')

methods=["constant box, merged 4 near boxes", "near neighbors"];

% methods=["constant box 2x", "constant box, merged 4 near boxes", "near neighbors"];
legendStrings=[];
for method = methods
    MethodName = sprintf('%s',method);
    path = ['E:\visualize\NewAnalysis\constant box Vs near neighbors\' MethodName '/']
    listL = readmatrix([path, 'list.txt'])
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

    DataIndex = 25; % acc
    % DataIndex = 30; % complexity
    % DataIndex = 31; % rnd
    if strcmp(method, "near neighbors")
        DataIndex = DataIndex + 1; 
    end
    
    Data = data(DataIndex,:);
    
    if DataIndex == 30
        if strcmp(method, "constant box, merged 4 near boxes")
            Data = data(30,:) - ((2*listL).^2)' .* (1-data(25,:));
        elseif strcmp(method, "constant box")
            Data = data(30,:) - ((listL).^2)' .* (1-data(25,:));
        elseif strcmp(method, "constant box 2x")
            Data = data(30,:) - ((2*listL).^2)' .* (1-data(25,:));
        end
    end
    plot(listL, Data, 'Marker', '.','MarkerSize', 18);
    hold on
    legendStrings = [legendStrings, method];
end

if DataIndex == 26
    % yline(0.238,'--','Metropolis')
    % ylim([0 0.4])
    ylim([0 0.3])
end
legendStrings=["constant boxes", "our work (only near neighbors box)"];

% legendStrings=["constant boxes with near neighbors in separate boxes", "constant boxes with near neighbors in one box", "our work (only near neighbors box)"];

% ylim([0 10])
yl = ylim;
xl = xlim;
% yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
% xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
set(gca,'fontsize',12) 
ylabel ('$P_{\mathrm{acc}}$'); name = 'acc';
% ylabel ('$N_\mathrm{PW}$'); name = 'pw';
% ylabel ('$N_\mathrm{RND}$'); name = 'rnd';
xlabel('$r_{\mathrm{box}}$')
% legendStrings = [legendStrings, "$P_{\mbox{Clock}}$"];;
% legendStrings = [legendStrings, "$P_{\mbox{OR-Clock}}$"];
% legendStrings = [""];
% hl = legend(legendStrings, 'Location','southwest','Interpreter','latex');
hl = legend(legendStrings, 'Location','northwest','Interpreter','latex');
set(hl,'Box','off')
% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex')
set(0,'defaulttextinterpreter','latex');

exportgraphics(gcf,[save_address '/', name '.png'],'Resolution',dpi')


