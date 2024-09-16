%% load data
clear all
clc
addpath(pwd);
addpath([pwd,'\matlab']);
addpath([pwd,'\matlab\mc']);
addpath([pwd,'\matlab\lattice2d']);
addpath([pwd,'\matlab\exportFig']);

path = ['E:/visualize/input/']
listL = readmatrix([path, 'list.txt'])
listN = listL.^2;
N = length(listN)

for i=1:N
    name = num2str(i);
    data_path = [path, name , '/'];

    % [status, msg, msgID] = mkdir(folder_address)

    data(:,:,i) = readmatrix([data_path 'MeanX.csv']);
    error(:,:,i) = readmatrix([data_path 'errorX.csv']);
end
% error = error / sqrt(10); % for results before 2022.08.25

date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
save_address = ['E:/visualize/output/', date_time_string];
[status, msg, msgID] = mkdir(save_address)

% myxLabel = '$k_{\mathrm{B}} T/\Lambda$';
myxLabel = '$T$';
% myxLabel = '$B\mu /\Lambda $';

dpi = 600;
fontSize = 14;

% orderParam = data(:,2,:);
orderParam = data(:,1,:);

addpath(pwd);
addpath([pwd,'\matlab\mc']);
addpath([pwd,'\matlab\lattice2d']);
%%---------------------------------------------------------------------------------
%% plot colors
mycolors = [];
for i=1:N
    mycolors = [mycolors;  rand(1,3)];
end
%%---------------------------------------------------------------------------------
%% Load locations & auto correlation data
% pathSammple0 = ['../visualize/input/', '0/']




%%---------------------------------------------------------------------------------
%% visualize ground state
% l = zeros(n ,1);
% A = zeros(n, 2);
% xr = zeros(n, 3);
% for j=1:3
%     for i=1:n
%         xr(i, j) = dipoles(3*i-3+j);
%     end
% end
% 
% for i = 1:n
% %     x = dipoles(1, i);
% %     y = dipoles(2, i);
% %     z = dipoles(3, i);
% %     
%     [A(i,2), A(i,1), l(i)] = cart2sph(xr(i, 1), xr(i, 2), xr(i, 3));
%     % second output is angle of elevation from xy plane
%     
%     A(i, 1) = pi/2 - A(i, 1); % convert to common spherical theta angle
%     
% end
% 
% l = l*0.15;
% 
% plotLatticeDipole(R, r, 1, A, l);
% axis off
% saveas(gcf, [folder_address '/', 'finalconfig.png'])
% saveas(gcf, [folder_address '/', 'finalconfig.fig'])
% clf('reset')
% 
% plotLattice(R, r, 1);
% saveas(gcf, [folder_address '/', 'positions.fig']);
% saveas(gcf, [folder_address '/', 'positions.png'])
% clf('reset')

%%---------------------------------------------------------------------------------
%% FSS Binder

Tc = 0.844;
reducedTemp = data(:,1,:) / Tc -1;
L = sqrt(listN);
for i=1:N
    Xtilde(:,i) = data(:,7,i);  %Binder
    x(:,i) = L(i)^(1) * reducedTemp(:,:,i);
end

clf('reset')
for i = 1:N
    plot(x(:,i), Xtilde(:,i) , '.','MarkerSize', 12);
%     plot(x(:,i), data(:,7,i), '.','MarkerSize', 10);
    hold on
end
% set(gca,'TickLabelInterpreter','latex');
% xlim([-5 5])
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));

% ylabel '$m L^{\beta / \nu}$'
% ylabel '$\chi L^{-\gamma / \nu}$'
ylabel '$U_4$'
xlabel('$t L^{1/ \nu}$')
legendStrings = "L = " + string(listL);
legend(legendStrings, 'Location','southeast','Interpreter','latex');
set(gca,'fontsize',fontSize) 
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
cmap = turbo(N+1);
ax.ColorOrder = cmap;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
exportgraphics(gcf,[save_address '/', 'Fss binder.png'],'Resolution',dpi)


%%---------------------------------------------------------------------------------
%% FSS chi

Tc = 0.844;
reducedTemp = data(:,1,:) / Tc -1;
L = sqrt(listN);
for i=1:N
    Xtilde(:,i) = L(i)^(-7/4) .* data(:,6,i);  
    x(:,i) = L(i)^(1) * reducedTemp(:,:,i);
end

clf('reset')
for i = 1:N
    plot(x(:,i), Xtilde(:,i) , '.','MarkerSize', 12);
%     plot(x(:,i), data(:,7,i), '.','MarkerSize', 10);
    hold on
end
% set(gca,'TickLabelInterpreter','latex');
% xlim([-5 5])
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
set(gca,'fontsize',fontSize) 
% ylabel '$m L^{\beta / \nu}$'
ylabel '$\chi L^{-\gamma / \nu}$'
% ylabel '$U_4$'
xlabel('$t L^{1/ \nu}$')
legendStrings = "L = " + string(listL);
legend(legendStrings, 'Location','northwest','Interpreter','latex');
set(gca,'fontsize',fontSize) 
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
cmap = turbo(N+1);
ax.ColorOrder = cmap;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
exportgraphics(gcf,[save_address '/', 'Fss chi.png'],'Resolution',dpi)


%%---------------------------------------------------------------------------------
%% FSS m

Tc = 0.844;
reducedTemp = data(:,1,:) / Tc -1;
L = sqrt(listN);
for i=1:N
    Xtilde(:,i) = L(i)^(1/8) .* data(:,5,i);  %Binder
    x(:,i) = L(i)^(1) * reducedTemp(:,:,i);
end

clf('reset')
for i = 1:N
    plot(x(:,i), Xtilde(:,i) , '.','MarkerSize', 12);
%     plot(x(:,i), data(:,7,i), '.','MarkerSize', 10);
    hold on
end
% set(gca,'TickLabelInterpreter','latex');
% xlim([-5 5])
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
set(gca,'fontsize',fontSize) 
ylabel '$m L^{\beta / \nu}$'
% ylabel '$\chi L^{-\gamma / \nu}$'
% ylabel '$U_4$'
xlabel('$t L^{1/ \nu}$')
legendStrings = "L = " + string(listL);
legend(legendStrings, 'Location','southwest','Interpreter','latex');
set(gca,'fontsize',fontSize) 
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
cmap = turbo(N+1);
ax.ColorOrder = cmap;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
exportgraphics(gcf,[save_address '/', 'Fss m.png'],'Resolution',dpi)

%%---------------------------------------------------------------------------------
%% E

clf('reset')
cmap = turbo(N+1);
for i = 1:N
    plot(orderParam(:,:,i), data(:,3,i), 'Marker', '.','MarkerSize', 15);
    hold on
end
% % whitebg; 
% set(gcf, 'Color','k');set(gcf, 'InvertHardCopy', 'off');

% set(gca, 'Color','k', 'XColor','w', 'YColor','w')
% set(gca,'TickLabelInterpreter','latex');
% mycolors = [];
% for i=1:N
%     mycolors = [mycolors;  rand(1,3)];
% end
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
set(gca,'fontsize',fontSize) 
ylabel '$E$'
% ylabel 'energy'
xlabel(myxLabel)
legendStrings = "L = " + string(listL);
hl = legend(legendStrings, 'Location','northwest','Interpreter','latex');
% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
ax.ColorOrder = cmap;
axis(ax, 'tight');
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
% set(gcf,'Color',[0 0 0]);
% saveas(gcf, [save_address '/', 'energy.png'])
exportgraphics(gcf,[save_address '/', 'energy.png'],'Resolution',dpi)
% exportgraphics(gcf,[save_address '/', 'energy.pdf'],'Resolution',dpi, ...
%     'BackgroundColor','none','ContentType','vector')
% print('[save_address '/', 'energy.pdf']','-dpng',,'-r600')

%%---------------------------------------------------------------------------------
%%  m
clf('reset')
set(gca,'XMinorTick','on','YMinorTick','on')
for i = 1:N
%     plot(orderParam(:,:,i), data(:,5,i), 'Marker', '.','MarkerSize', 15);
    errorbar(orderParam(:,:,i), data(:,5,i), error(:,5,i), 'Marker', '.','MarkerSize', 15);
    hold on
end
% whitebg; set(gcf, 'Color','k');set(gcf, 'InvertHardCopy', 'off')
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
xlabel(myxLabel)
ylabel '$m$'
% ylabel 'magnetization'
legendStrings = "L = " + string(listL)
legend(legendStrings, 'Location','southwest','Interpreter','latex');
% legend(legendStrings, 'Location','southwest','Interpreter','latex','fontsize',8);
set(gca,'fontsize',fontSize) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
ax.ColorOrder = cmap;
axis(ax, 'tight') 
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
exportgraphics(gcf,[save_address '/', 'm.png'],'Resolution',dpi)
% saveas(gcf, [save_address '/', 'm.png'])

%%---------------------------------------------------------------------------------
%% C
clf('reset')
set(gca,'XMinorTick','on','YMinorTick','on')
for i = 1:N
    errorbar(orderParam(:,:,i), data(:,4,i), error(:,4,i), 'Marker', '.','MarkerSize', 15);
%     plot(orderParam(:,:,i), data(:,4,i), 'Marker', '.','MarkerSize', 15);
    hold on
end
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
xlabel(myxLabel)
ylabel '$c$'
% ylabel 'heat capaity'
legendStrings = "L = " + string(listL)
% legend(legendStrings,'Interpreter','latex')
legend(legendStrings, 'Location','northwest','Interpreter','latex');

set(gca,'fontsize',fontSize) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
ax.ColorOrder = cmap;
axis(ax, 'tight');

% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontsize',18)  
% axes('Position',[.65 .65 .25 .25]) % begin inset
% box on
% scatter(T_E_m_c_X(:,1), T_E_m_c_X(:,2), 15, 'filled');
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontsize',fontSize) 
% ylabel '$E$'
% xlabel '$T$'
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
exportgraphics(gcf,[save_address '/', 'C.png'],'Resolution',dpi)
% saveas(gcf, [save_address '/', 'C.png']) % end inset

%%---------------------------------------------------------------------------------
%% chi
clf('reset')
set(gca,'XMinorTick','on','YMinorTick','on')
for i = 1:N
    errorbar(orderParam(:,:,i), data(:,6,i), error(:,6,i), 'Marker', '.','MarkerSize', 15);
%     plot(orderParam(:,:,i), data(:,6,i), 'Marker', '.','MarkerSize', 15);
    hold on
end
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
xlabel(myxLabel)
ylabel '$\chi$'
% ylabel 'susceptibility'
legendStrings = "L = " + string(listL)
legend(legendStrings,'Interpreter','latex')
set(gca,'fontsize',fontSize) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
ax.ColorOrder = cmap;
axis(ax, 'tight');
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
exportgraphics(gcf,[save_address '/', 'X.png'],'Resolution',dpi)
% saveas(gcf, [save_address '/', 'X.png'])
%%---------------------------------------------------------------------------------
%% Binder
% clf('reset')
f = figure;
f.Position = [100 100 300 200];
xRange = linspace(data(end,1,1), data(1,1,1), 50); 
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', ...
           'x', '+', '*', 'o', 'x', 's', 'd', '^', 'v', '>','<', 'p', 'h'};

LIndicesPlot = [3,5,6,8,10];
cmap = linspecer(length(LIndicesPlot)+1);
mi = 1;
for i = LIndicesPlot
%     SplineVal = spline(orderParam(:,:,i), data(:,6,i),xRange);
%     errorbar(orderParam(:,:,i), data(:,6,i), error(:,6,i), '.', xRange, SplineVal);
%     f = fit(orderParam(:,:,i),data(:,6,i),'smoothingspline');
    % DataToPlot = 1 - data(:,7,i) ./ 3;
    % ErrorPlot = error(:,7,i) ./ 3;
    DataToPlot = data(:,7,i);
    ErrorPlot = error(:,7,i);
    % errorbar(orderParam(:,:,i), DataToPlot, ErrorPlot, '--.','MarkerSize', 15);
    h = errorbar(orderParam(:,:,i), DataToPlot, ErrorPlot, ...
             'Color',cmap(mi,:), 'Marker', markers{mi}','MarkerSize', 5);
    h.MarkerFaceColor = h.Color;
%     plot(orderParam(:,:,i), data(:,7,i), '--.','MarkerSize', 15);
    hold on
%     plot(f); 
    mi = mi + 1;
end
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/4) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/4) :xl(2));

% axis tight
xlabel(myxLabel)
ylabel '$U_4$'
legendStrings = "L = " + string(listL(LIndicesPlot))
lgd = legend(legendStrings, 'Location','southwest','Interpreter','latex');
lgd.FontSize = 11;
set(lgd,'Box','off')
set(gca,'fontsize',13) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

% axes('Position',[0.2 .65 .25 .25]) % begin inset
% box on
% for i = 1:N
%     plot(orderParam(:,:,i), data(:,7,i), '--.','MarkerSize', 9);
%     hold on
% end
% xlim([0.840 0.850])
% ylim([1.06 1.14])
% % ylim([1 1.2])
% % xticks(0.64:0.02:0.68)
% % xlabel '$T$'
% % ylabel '$Binder parameter$'
% ax = gca;
% ax.TickLength(1) = 0.04;
% ax.ColorOrder = cmap;
% 
% % ax.YAxis.MinorTickValues = ax.YAxis.Limits(1):0.25:ax.YAxis.Limits(2);
% % ax.XAxis.MinorTickValues = ax.XAxis.Limits(1):0.5:ax.XAxis.Limits(2);
% set(gca,'XMinorTick','on','YMinorTick','on')
% % set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontsize',9)  % end inset
% set(gca,'TickLabelInterpreter','latex');
% set(0,'defaulttextinterpreter','latex');
exportgraphics(gcf,[save_address '/', 'Binder.png'],'Resolution',dpi, ...
    'BackgroundColor','none')
% % saveas(gcf, [save_address '/', 'Binder.png'])

%%---------------------------------------------------------------------------------
%% chi max vs N, power fit
clf('reset')
% set(gca,'XMinorTick','on','YMinorTick','on')
invsize = 1 ./ listN;
sizes = sqrt(listN);
for i = 1:N
%     TT = 0;
%     chidiff = 0;
    [maxchi, maxidx] = max(data(:,6,i));
    maxidxL(i) = maxchi;
%     TT(1:maxidx-1) = data(1:maxidx-1, 1, i);
%     TT(maxidx+1:1) = data(maxidx+1:1, 1, i);
    
%     chidiff(1:maxidx-1) = maxchi - data(1:maxidx-1, 5, i);
%     chidiff(maxidx+1:1) = maxchi - data(maxidx+1:1, 5, i);
%     loglog(TT, chidiff,'Marker', '.','MarkerSize', 15);
%     loglog(TT, chidiff);
%     hold on
end
sizes = sizes';
maxidxL = maxidxL';
% f = fit(sizes,maxidxL,'exp1')
plot(sizes, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
hold on
% % logsizes = log(sizes);
% logmaxidxL = log(maxidxL);
% plot(logsizes, logmaxidxL, 's', 'Marker', '.','MarkerSize', 15);
% 
% plot(sizes, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
% p = polyfit(log(sizes),maxidxL,1);
% lineStartEnd = [1, max(sizes)+10]
% z = polyval(p,log(sizes)); 
% hold on;
% plot(sizes,z);

xFit = linspace(min(sizes) - min(sizes)/10, max(sizes) + max(sizes)/10, 500);
p = polyfit(log(sizes),log(maxidxL),1)
yFit = polyval(p, log(xFit));
plot(xFit, exp(yFit) ,'r' , 'LineWidth', 1 );

xl = xlim;
yl = ylim;
xt = 0.1 * (xl(2)-xl(1)) + xl(1);
yt = 0.2 * (yl(2)-yl(1)) + yl(1);
caption = sprintf('$\\chi \\sim  L^{%.2f}$', p(1));
h = text(xt, yt, caption, 'FontSize', fontSize, 'Color', 'r', 'FontWeight', 'bold');
set(h,'Rotation',30);

% axis padded
xlabel '$N$'
ylabel '$\chi_{\mathrm{max}}$'
ylim([1 (max(maxidxL)+ max(maxidxL)/3)])
set(gca, 'YScale', 'log', 'XScale', 'log')   %using plot
% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
% set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'fontsize',12) 
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
exportgraphics(gcf,[save_address '/', 'XvsN power fit.png'],'Resolution',dpi, ...
    'BackgroundColor','none')

%%---------------------------------------------------------------------------------
%% chi max vs N, log fit
clf('reset')
% set(gca,'XMinorTick','on','YMinorTick','on')
invsize = 1 ./ listN;
sizes = sqrt(listN);
for i = 1:N
%     TT = 0;
%     chidiff = 0;
    [maxchi, maxidx] = max(data(:,6,i));
    maxidxL(i) = maxchi;
%     TT(1:maxidx-1) = data(1:maxidx-1, 1, i);
%     TT(maxidx+1:1) = data(maxidx+1:1, 1, i);
    
%     chidiff(1:maxidx-1) = maxchi - data(1:maxidx-1, 5, i);
%     chidiff(maxidx+1:1) = maxchi - data(maxidx+1:1, 5, i);
%     loglog(TT, chidiff,'Marker', '.','MarkerSize', 15);
%     loglog(TT, chidiff);
%     hold on
end
sizes = sizes';
maxidxL = maxidxL';
% f = fit(sizes,maxidxL,'exp1')
plot(sizes, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
hold on
% % logsizes = log(sizes);
% logmaxidxL = log(maxidxL);
% plot(logsizes, logmaxidxL, 's', 'Marker', '.','MarkerSize', 15);
% 
% plot(sizes, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
% p = polyfit(log(sizes),maxidxL,1);
% lineStartEnd = [1, max(sizes)+10]
% z = polyval(p,log(sizes)); 
% hold on;
% plot(sizes,z);

xFit = linspace(min(sizes) - min(sizes)/10, max(sizes)+ max(sizes)/10, 500);
p = polyfit(log(sizes),maxidxL,1)
yFit = polyval(p, log(xFit));
plot(xFit, yFit ,'r' , 'LineWidth', 1 );

xl = xlim;
yl = ylim;
xt = 0.1 * (xl(2)-xl(1)) + xl(1);
yt = 0.5 * (yl(2)-yl(1)) + yl(1);
caption = sprintf('$\\chi \\sim  \\log{L}$');
h = text(xt, yt, caption, 'FontSize', fontSize, 'Color', 'r', 'FontWeight', 'bold');
set(h,'Rotation',42);

xlabel '$N$'
ylabel 'Maximum of susceptibility'
% ylim([-10 sizes(end)+10])
% set(gca, 'YScale', 'log', 'XScale', 'log')   %using plot
set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
% set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'fontsize',12) 
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'XvsN Log fit.png'])

%%---------------------------------------------------------------------------------
%% c max vs N, log fit
clf('reset')
% set(gca,'XMinorTick','on','YMinorTick','on')
invsize = 1 ./ listN;
sizes = sqrt(listN);
for i = 1:N
%     TT = 0;
%     chidiff = 0;
    [maxchi, maxidx] = max(data(:,4,i));
    maxidxL(i) = maxchi;
%     TT(1:maxidx-1) = data(1:maxidx-1, 1, i);
%     TT(maxidx+1:1) = data(maxidx+1:1, 1, i);
    
%     chidiff(1:maxidx-1) = maxchi - data(1:maxidx-1, 5, i);
%     chidiff(maxidx+1:1) = maxchi - data(maxidx+1:1, 5, i);
%     loglog(TT, chidiff,'Marker', '.','MarkerSize', 15);
%     loglog(TT, chidiff);
%     hold on
end
sizes = sizes';
maxidxL = maxidxL';
% f = fit(sizes,maxidxL,'exp1')
plot(sizes, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
hold on
% % logsizes = log(sizes);
% logmaxidxL = log(maxidxL);
% plot(logsizes, logmaxidxL, 's', 'Marker', '.','MarkerSize', 15);
% 
% plot(sizes, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
% p = polyfit(log(sizes),maxidxL,1);
% lineStartEnd = [1, max(sizes)+10]
% z = polyval(p,log(sizes)); 
% hold on;
% plot(sizes,z);

xFit = linspace(min(sizes) - min(sizes)/10, max(sizes)+ max(sizes)/10, 500);
p = polyfit(log(sizes),maxidxL,1)
yFit = polyval(p, log(xFit));
plot(xFit, yFit ,'r' , 'LineWidth', 1 );

xl = xlim;
yl = ylim;
xt = 0.2 * (xl(2)-xl(1)) + xl(1);
yt = 0.4 * (yl(2)-yl(1)) + yl(1);
caption = sprintf('$c \\sim  \\log{L}$');
h = text(xt, yt, caption, 'FontSize', fontSize, 'Color', 'r', 'FontWeight', 'bold');
set(h,'Rotation',42);

xlabel '$N$'
ylabel 'Maximum of heat capacity'
% ylim([-10 sizes(end)+10])
% set(gca, 'YScale', 'log', 'XScale', 'log')   %using plot
set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
% set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'fontsize',12) 
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'CvsN Log fit.png'])

%%---------------------------------------------------------------------------------
%% c max vs N, power fit 
clf('reset')
% set(gca,'XMinorTick','on','YMinorTick','on')
invsize = 1 ./ listN;
sizes = sqrt(listN);
for i = 1:N
%     TT = 0;
%     chidiff = 0;
    [maxchi, maxidx] = max(data(:,4,i));
    maxidxL(i) = maxchi;
%     TT(1:maxidx-1) = data(1:maxidx-1, 1, i);
%     TT(maxidx+1:1) = data(maxidx+1:1, 1, i);
    
%     chidiff(1:maxidx-1) = maxchi - data(1:maxidx-1, 5, i);
%     chidiff(maxidx+1:1) = maxchi - data(maxidx+1:1, 5, i);
%     loglog(TT, chidiff,'Marker', '.','MarkerSize', 15);
%     loglog(TT, chidiff);
%     hold on
end
sizes = sizes';
maxidxL = maxidxL';
% f = fit(sizes,maxidxL,'exp1')
plot(sizes, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
hold on
% % logsizes = log(sizes);
% logmaxidxL = log(maxidxL);
% plot(logsizes, logmaxidxL, 's', 'Marker', '.','MarkerSize', 15);
% 
% plot(sizes, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
% p = polyfit(log(sizes),maxidxL,1);
% lineStartEnd = [1, max(sizes)+10]
% z = polyval(p,log(sizes)); 
% hold on;
% plot(sizes,z);

xFit = linspace(min(sizes) - min(sizes)/10, max(sizes) + max(sizes)/10, 500);
p = polyfit(log(sizes),log(maxidxL),1)
yFit = polyval(p, log(xFit));
plot(xFit, exp(yFit) ,'r' , 'LineWidth', 1 );

xl = xlim;
yl = ylim;
xt = 0.1 * (xl(2)-xl(1)) + xl(1);
yt = -0.5 * (yl(2)-yl(1)) + yl(1);
caption = sprintf('$c \\sim  L^{%.2f}$', p(1));
h = text(xt, yt, caption, 'FontSize', fontSize, 'Color', 'r', 'FontWeight', 'bold');
set(h,'Rotation',0);

% axis padded
xlabel '$N$'
ylabel 'Maximum of heat capacity'
ylim([1 (max(maxidxL)+ max(maxidxL)/3)])
set(gca, 'YScale', 'log', 'XScale', 'log')   %using plot
% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
% set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'fontsize',12) 
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'CvsN power fit.png'])

%%---------------------------------------------------------------------------------
%% C max vs N, log log fit
clf('reset')
% set(gca,'XMinorTick','on','YMinorTick','on')
invsize = 1 ./ listN;
sizes = listN;
for i = 1:N
%     TT = 0;
%     chidiff = 0;
    [maxchi, maxidx] = max(data(:,4,i));
    maxidxL(i) = maxchi;
%     TT(1:maxidx-1) = data(1:maxidx-1, 1, i);
%     TT(maxidx+1:1) = data(maxidx+1:1, 1, i);
    
%     chidiff(1:maxidx-1) = maxchi - data(1:maxidx-1, 5, i);
%     chidiff(maxidx+1:1) = maxchi - data(maxidx+1:1, 5, i);
%     loglog(TT, chidiff,'Marker', '.','MarkerSize', 15);
%     loglog(TT, chidiff);
%     hold on
end
sizes = sizes';
maxidxL = maxidxL';
% f = fit(sizes,maxidxL,'exp1')
plot(sizes, log(log(maxidxL)), 's', 'Marker', '.','MarkerSize', 15);
hold on
% % logsizes = log(sizes);
% logmaxidxL = log(maxidxL);
% plot(logsizes, logmaxidxL, 's', 'Marker', '.','MarkerSize', 15);
% 
% plot(sizes, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
% p = polyfit(log(sizes),maxidxL,1);
% lineStartEnd = [1, max(sizes)+10]
% z = polyval(p,log(sizes)); 
% hold on;
% plot(sizes,z);

% xl = xlim;
% yl = ylim;
% xt = 0.2 * (xl(2)-xl(1)) + xl(1);
% yt = 0.5 * (yl(2)-yl(1)) + yl(1);
% caption = sprintf('$c \\sim  \\log{N}$');
% h = text(xt, yt, caption, 'FontSize', fontSize, 'Color', 'r', 'FontWeight', 'bold');
% set(h,'Rotation',42);

xlabel '$N$'
ylabel '$\log(\log(C_{\mathrm{max}}))$'
% ylim([-10 sizes(end)+10])
% set(gca, 'YScale', 'log', 'XScale', 'log')   %using plot
% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
% set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'fontsize',12) 
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'CvsN Log log fit.png'])


%%---------------------------------------------------------------------------------
%%  m plane
clf('reset')
set(gca,'XMinorTick','on','YMinorTick','on')
for i = 1:N
    plot(orderParam(:,:,i), data(:,8,i), 'Marker', '.','MarkerSize', 15);
    hold on
end
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
xlabel(myxLabel)
ylabel '$|\textbf{m}_{\mathrm{plane}}|$'
legendStrings = "N = " + string(listN)
legend(legendStrings, 'Location','southwest','Interpreter','latex')
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'm plane.png'])
%%---------------------------------------------------------------------------------
%% chi plane
clf('reset')
set(gca,'XMinorTick','on','YMinorTick','on')
for i = 1:N
%     errorbar(orderParam(:,:,i), data(:,5,i), error(:,5,i), 'Marker', '.','MarkerSize', 15);
    plot(orderParam(:,:,i), data(:,9,i), 'Marker', '.','MarkerSize', 15);
    hold on
end
xlabel(myxLabel)
ylabel '$\chi_{\mathrm{plane}}$'
legendStrings = "N = " + string(listN)
legend(legendStrings,'Interpreter','latex')
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'X plane.png'])
%%---------------------------------------------------------------------------------
%% Binder plane
clf('reset')
fig = figure();
ax = axes(fig);
set(gca,'XMinorTick','on','YMinorTick','on')
xRange = linspace(data(end,1,1), data(1,1,1), 50); 
for i = 3:N
%     SplineVal = spline(orderParam(:,:,i), data(:,6,i),xRange);
%     errorbar(orderParam(:,:,i), data(:,6,i), error(:,6,i), '.', xRange, SplineVal);
%     f = fit(orderParam(:,:,i),data(:,6,i),'smoothingspline');
%     errorbar(orderParam(:,:,i), data(:,6,i), error(:,6,i), '--.','MarkerSize', 15);
    plot(orderParam(:,:,i), data(:,10,i), '--.','MarkerSize', 15);
    hold on
%     plot(f);    
end
xlabel(myxLabel)
ylabel '$U_{\mathrm{plane}}$'
legendStrings = "N = " + string(listN)
lgd = legend(legendStrings, 'Location','southeast','Interpreter','latex')
lgd.FontSize = 9;
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

axes('Position',[0.2 .65 .25 .25]) % begin inset
box on
for i = 1:N
    plot(orderParam(:,:,i), data(:,9,i), '--.','MarkerSize', 9);
    hold on
end
xlim([0.64 0.7])
ylim([1.04 1.2])
% ylim([1.01 1.025])
% xticks(0.64:0.02:0.68)
% xlabel '$T$'
% ylabel '$Binder parameter$'
ax = gca;
ax.TickLength(1) = 0.04;
% ax.YAxis.MinorTickValues = ax.YAxis.Limits(1):0.25:ax.YAxis.Limits(2);
% ax.XAxis.MinorTickValues = ax.XAxis.Limits(1):0.5:ax.XAxis.Limits(2);
set(gca,'XMinorTick','on','YMinorTick','on')
% set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',9)  % end inset
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'Binder plane.png'])

%%---------------------------------------------------------------------------------
%% chi plane max vs N, power fit
clf('reset')
% set(gca,'XMinorTick','on','YMinorTick','on')
invsize = 1 ./ listN;
sizes = listN;
for i = 1:N
%     TT = 0;
%     chidiff = 0;
    [maxchi, maxidx] = max(data(:,8,i));
    maxidxL(i) = maxchi;
%     TT(1:maxidx-1) = data(1:maxidx-1, 1, i);
%     TT(maxidx+1:1) = data(maxidx+1:1, 1, i);
    
%     chidiff(1:maxidx-1) = maxchi - data(1:maxidx-1, 5, i);
%     chidiff(maxidx+1:1) = maxchi - data(maxidx+1:1, 5, i);
%     loglog(TT, chidiff,'Marker', '.','MarkerSize', 15);
%     loglog(TT, chidiff);
%     hold on
end
sizes = sizes';
maxidxL = maxidxL';
% f = fit(sizes,maxidxL,'exp1')
plot(sizes, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
hold on
% % logsizes = log(sizes);
% logmaxidxL = log(maxidxL);
% plot(logsizes, logmaxidxL, 's', 'Marker', '.','MarkerSize', 15);
% 
% plot(sizes, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
% p = polyfit(log(sizes),maxidxL,1);
% lineStartEnd = [1, max(sizes)+10]
% z = polyval(p,log(sizes)); 
% hold on;
% plot(sizes,z);

xFit = linspace(min(sizes) - min(sizes)/10, max(sizes) + max(sizes)/10, 500);
p = polyfit(log(sizes),log(maxidxL),1)
yFit = polyval(p, log(xFit));
plot(xFit, exp(yFit) ,'r' , 'LineWidth', 1 );

xl = xlim;
yl = ylim;
xt = 0.1 * (xl(2)-xl(1)) + xl(1);
yt = 0.2 * (yl(2)-yl(1)) + yl(1);
caption = sprintf('$\\chi \\sim  N^{%.2f}$', p(1));
h = text(xt, yt, caption, 'FontSize', fontSize, 'Color', 'r', 'FontWeight', 'bold');
set(h,'Rotation',30);

% axis padded
xlabel '$N$'
ylabel '$\chi_{\mathrm{plane, max}}$'
ylim([1 (max(maxidxL)+ max(maxidxL)/3)])
set(gca, 'YScale', 'log', 'XScale', 'log')   %using plot
% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
% set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
saveas(gcf, [save_address '/', 'Xplane vs N power fit.png'])

%%---------------------------------------------------------------------------------
%%  mz
clf('reset')
set(gca,'XMinorTick','on','YMinorTick','on')
for i = 1:N
    plot(orderParam(:,:,i), data(:,13,i), 'Marker', '.','MarkerSize', 15);
    hold on
end
yl = ylim;
xl = xlim;
yticks(yl(1): ((yl(2)- yl(1))/5) :yl(2));
xticks(xl(1): ((xl(2)- xl(1))/5) :xl(2));
xlabel(myxLabel)
ylabel '${m}_{\mathrm{z}}$'
legendStrings = "N = " + string(listN)
legend(legendStrings, 'Location','northwest','Interpreter','latex')
set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
saveas(gcf, [save_address '/', 'mz.png'])

%%---------------------------------------------------------------------------------
%% read cursor loction
 button = 1;
 while sum(button) <=3   % read ginputs until a mouse right-button occurs
   [x,y,button] = ginput(3);
 end
%%---------------------------------------------------------------------------------
%% loglog chi (not usable)
clf('reset')
% set(gca,'XMinorTick','on','YMinorTick','on')
for i = 1:N
    TT = 0;
    chidiff = 0;
    [maxchi, maxidx] = max(data(:,5,i))
    TT(1:maxidx-1) = data(1:maxidx-1, 1, i);
    TT(maxidx+1:1) = data(maxidx+1:1, 1, i);
    
    chidiff(1:maxidx-1) = maxchi - data(1:maxidx-1, 5, i);
    chidiff(maxidx+1:1) = maxchi - data(maxidx+1:1, 5, i);
    loglog(TT, chidiff,'Marker', '.','MarkerSize', 15);
%     loglog(TT, chidiff);
    hold on
end
xlabel(myxLabel)
ylabel '$\chi - \chi_{\mathrm{max}}$'
legendStrings = "N = " + string(listN)
legend(legendStrings)
% set(gca, 'YScale', 'log', 'XScale', 'log')   %using plot
% set(gca,'fontsize',12) 
% set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
saveas(gcf, [save_address '/', 'Xloglog.png'])

%%---------------------------------------------------------------------------------
%% loglog C (not usable)
clf('reset')
% set(gca,'XMinorTick','on','YMinorTick','on')
for i = 1:N
    TT = zeros(1,1);
    chidiff = zeros(1,1);
    [maxchi, maxidx] = max(data(:,4,i))
    TT(1:maxidx-1) = data(1:maxidx-1, 1, i);
    TT(maxidx+1:1) = data(maxidx+1:1, 1, i);
    
    chidiff(1:maxidx-1) = maxchi - data(1:maxidx-1, 4, i);
    chidiff(maxidx+1:1) = maxchi - data(maxidx+1:1, 4, i);
    loglog(TT, chidiff,'Marker', '.','MarkerSize', 15);
    hold on
end
xlabel(myxLabel)
ylabel '$C - C_{\mathrm{max}}$'
legendStrings = "N = " + string(listN)
legend(legendStrings)
set(gca,'fontsize',12) 
% set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
saveas(gcf, [save_address '/', 'Cloglog.png'])

%%---------------------------------------------------------------------------------
%% temp vs 1/L, fit 

clf('reset')
% set(gca,'XMinorTick','on','YMinorTick','on')
invsize = 1 ./ listN;
% sizes = 1./ [8,16,24,32,48,64] ;
sizes = invsize ;
jexcl = 0; %% number of excluded from fit
TcNbyChi = zeros(1,N-jexcl);
TcNbyC = zeros(1,N-jexcl);
for i = (1+jexcl):N
%     TT = 0;
%     chidiff = 0;
    [~, maxidxChi] = max(data(:,5,i));
    [~, maxidxC] = max(data(:,4,i));
    TcNbyChi(i-jexcl) = data(maxidxChi, 1, i);
    TcNbyC(i-jexcl)   = data(maxidxC, 1, i);
%     TT(1:maxidx-1) = data(1:maxidx-1, 1, i);
%     TT(maxidx+1:1) = data(maxidx+1:1, 1, i);
    
%     chidiff(1:maxidx-1) = maxchi - data(1:maxidx-1, 5, i);
%     chidiff(maxidx+1:1) = maxchi - data(maxidx+1:1, 5, i);
%     loglog(TT, chidiff,'Marker', '.','MarkerSize', 15);
%     loglog(TT, chidiff);
%     hold on
end
sizes = sizes';
TcNbyChi = TcNbyChi';
TcNbyC = TcNbyC';
% f = fit(sizes,maxidxL,'poly1')
mdlchi = fitlm(sizes,TcNbyChi);
betahatchi =mdlchi.Coefficients.Estimate
mdlc = fitlm(sizes,TcNbyC);
betahatc =mdlc.Coefficients.Estimate
% anova(mdl,'summary')
plot(sizes, TcNbyChi, 's' ,'MarkerEdgeColor','blue',...
    'MarkerFaceColor', 'blue', 'MarkerSize', 9 )
hold on
xFit = linspace(0, max(sizes)+ max(sizes)/10, 500);
p = polyfit(sizes, TcNbyChi, 1)
yFit = polyval(p, xFit);
plot(xFit, yFit ,'b' , 'LineWidth', 1 );

plot(sizes, TcNbyC, 'r.' ,'MarkerSize', 20)
hold on
xFit = linspace(0, max(sizes)+ max(sizes)/10, 500);
p = polyfit(sizes, TcNbyC, 1)
yFit = polyval(p, xFit);
plot(xFit, yFit ,'r', 'LineWidth', 1);
lgd = legend('$\chi_{max}$ locations', ' fit of $\chi$', ...
            '$c_{max}$ locations', 'fit of $c$', ...
            'Location','southeast', 'Interpreter','latex');
lgd.FontSize = 10;
% plt1 = plot(mdlchi);

% plt2 = plot(mdlc);
% plot(sizes, maxidxL,'Marker', '.','MarkerSize', 15)
% plt1.Matker = '.';
% logsizes = log(sizes);
% logmaxidxL = log(maxidxL);
% plot(logsizes, logmaxidxL, 's', 'Marker', '.','MarkerSize', 15);
% 
% plot(invsize, maxidxL, 's', 'Marker', '.','MarkerSize', 15);
% p = polyfit(log(invsize),log(maxidxL),1);
% lineStartEnd = [1, max(invsize)+10]
% z = polyval(p,log(invsize)); 
% hold on;
% plot(invsize,exp(z));
% polyval(p,log(1/1000))

xlabel '$N^{-1}$'
ylabel '$Tc(N)$'
% ylim([0.66 max(maxidxL)+0.01])
% xlim([0.01 max(invsize)+0.01])
% set(gca, 'YScale', 'log')   %using plot
% set(gca, 'YScale', 'log', 'XScale', 'log')   %using plot
% set(gca,'fontsize',12) 
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.02;
saveas(gcf, [save_address '/', 'TcvsL.png'])



%%   chi plot with ticks
%%---------------------------------------------------------------------------------
fig1 = plot(T_E_m_c_X(:,1), T_E_m_c_X(:,5), '.', 'MarkerSize', 20)
xlabel '$T$'
ylabel '$\chi$'
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.03;
ax.YAxis.MinorTickValues = ax.YAxis.Limits(1):0.25:ax.YAxis.Limits(2);
% ax.XAxis.MinorTickValues = ax.XAxis.Limits(1):0.25:ax.XAxis.Limits(2);
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',18)  
% axes('Position',[.53 .55 .35 .35]) % begin inset
% box on
% fig1 = plot(T_E_m_c_X(:,1), T_E_m_c_X(:,3), '.', 'MarkerSize', 15);
% xlabel '$T$'
% ylabel '$m$'
% ax = gca;
% ax.TickLength(1) = 0.04;
% % ax.YAxis.MinorTickValues = ax.YAxis.Limits(1):0.25:ax.YAxis.Limits(2);
% % ax.XAxis.MinorTickValues = ax.XAxis.Limits(1):0.5:ax.XAxis.Limits(2);
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontsize',fontSize)  % end inset
saveas(gcf, [folder_address '/', 'X.png'])
close()