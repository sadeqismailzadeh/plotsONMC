%% save address
clear all
clc
addpath(pwd);
addpath([pwd,'\matlab\mc']);
addpath([pwd,'\matlab\lattice2d']);

path = ['E:/visualize/inputHistogram/histogram/'];
date_time_string = datestr(now,'yyyy.mm.dd HH-MM-SS');
save_address = ['E:/visualize/output/', date_time_string , '_r2plus1r'];
[status, msg, msgID] = mkdir(save_address)

dpi = 600;

%%---------------------------------------------------------------------------------
%% plot


fig = figure;
fig.Position = [100 100 390 220];

legendStrings=[];

% Define a range of values for r
r = linspace(2, 30, 1000); % Adjust the range as needed, avoiding r = 0

% Calculate the function values
y = 0.23 * (1*pi*r.^2 + 16*pi*(1/0.849)*20*1./r);

% Plot the function
plot(r, y , 'b', 'LineWidth', 0.8)
legendStrings = [legendStrings, string(['$T = ' num2str(0.849), '$'])];
hold on
% Find the minimum value of the function
[min_y, min_index] = min(y);
min_r = r(min_index);
% Mark the minimum point on the plot
plot(min_r, min_y, 'bo', 'MarkerFaceColor', 'blue', 'HandleVisibility', 'off') % Mark the minimum point with a red circle



% Calculate the function values
y = 0.23 *  (1*pi*r.^2 + 16*pi*(1/0.2)*20*1./r);
plot(r, y, '--r', 'LineWidth', 0.8)
legendStrings = [legendStrings, string(['$T = ' num2str(0.2), '$'])];
% Find the minimum value of the function
[min_y, min_index] = min(y);
min_r = r(min_index);
% Mark the minimum point on the plot
plot(min_r, min_y, 'ro', 'MarkerFaceColor', 'red', 'HandleVisibility', 'off') % Mark the minimum point with a red circle


myFontSize =13;
set(gca,'fontsize',myFontSize) 

ylabel ('$t_\mathrm{spin}$', 'FontSize', myFontSize)   % show y label
xlabel('$r_\mathrm{box}$', 'FontSize', myFontSize) % show x label
leg1 = legend(legendStrings, 'Location','southeast','Interpreter','latex', 'FontSize', 12);
set(leg1,'Box','off')

% pbaspect([2 1 1])
% set(hl, 'TextColor','w')
set(gca,'XMinorTick','on','YMinorTick','on');
ax = gca;
ax.TickLength(1) = 0.02;
set(gca,'TickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');

% set(gca, 'XScale', 'log')   %using plot
% set(gca, 'YScale', 'log')   %using plot
exportgraphics(gcf,[save_address '/',  'r2plus1r' '.png'],'Resolution',dpi)
exportgraphics(gcf,[save_address '/',  'r2plus1r' '.eps'])
