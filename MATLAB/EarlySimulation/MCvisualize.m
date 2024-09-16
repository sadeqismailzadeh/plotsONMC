%% load data

% name = 'triangular'
% folder_address = ['../results/', name];
% [status, msg, msgID] = mkdir(folder_address)

addpath(pwd);
addpath([pwd,'\mc']);
addpath([pwd,'\lattice2d']);
%% visualize ground state

l = zeros(n ,1);
A = zeros(n, 2);
xr = zeros(n, 3);
for j=1:3
    for i=1:n
        xr(i, j) = dipoles(3*i-3+j);
    end
end

for i = 1:n
%     x = dipoles(1, i);
%     y = dipoles(2, i);
%     z = dipoles(3, i);
%     
    [A(i,2), A(i,1), l(i)] = cart2sph(xr(i, 1), xr(i, 2), xr(i, 3));
    % second output is angle of elevation from xy plane
    
    A(i, 1) = pi/2 - A(i, 1); % convert to common spherical theta angle
    
end

l = l*0.15;

plotLatticeDipole(R, r, 1, A, l);
axis off
saveas(gcf, [folder_address '/', 'finalconfig.png'])
saveas(gcf, [folder_address '/', 'finalconfig.fig'])
clf('reset')

plotLattice(R, r, 1);
saveas(gcf, [folder_address '/', 'positions.fig']);
saveas(gcf, [folder_address '/', 'positions.png'])
clf('reset')
%% visualize statistics
% set(gca,'TickLabelInterpreter','latex');

set(0,'defaulttextinterpreter','latex')

 
% scatter(T_E_C_op(:,1), T_E_C_op(:,2), 25, 'filled');
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontsize',18) 
% ylabel '$E$'
% xlabel '$T$' 
% saveas(gcf, [folder_address '/', 'energy.png'])
% clf('reset')
% 
% fig1 = scatter(T_E_C_op(:,1), T_E_C_op(:,4), 25, 'filled');
% xlabel '$T$'
% ylabel '$m$'
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontsize',18)  
% saveas(gcf, [folder_address '/', 'orderparam.png'])
% clf('reset')
% 
% fig1 = scatter(T_E_C_op(:,1), T_E_C_op(:,3), 25, 'filled');
% xlabel '$T$'
% ylabel '$C$'
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontsize',18)  
% axes('Position',[.65 .65 .25 .25]) % begin inset
% box on
% scatter(T_E_C_op(:,1), T_E_C_op(:,2), 15, 'filled');
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontsize',14) 
% ylabel '$E$'
% xlabel '$T$' 
% saveas(gcf, [folder_address '/', 'C.png']) % end inset
% clf('reset')

%-------------- plot X ----------------%
fig1 = plot(T_E_C_op(:,1), T_E_C_op(:,5), '.', 'MarkerSize', 20)
xlabel '$T$'
ylabel '$\chi$'
set(gca,'XMinorTick','on','YMinorTick','on')
ax = gca;
ax.TickLength(1) = 0.03;
ax.YAxis.MinorTickValues = ax.YAxis.Limits(1):0.25:ax.YAxis.Limits(2);
% ax.XAxis.MinorTickValues = ax.XAxis.Limits(1):0.25:ax.XAxis.Limits(2);
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',18)  
axes('Position',[.53 .55 .35 .35]) % begin inset
box on
fig1 = plot(T_E_C_op(:,1), T_E_C_op(:,4), '.', 'MarkerSize', 15);
xlabel '$T$'
ylabel '$m$'
ax = gca;
ax.TickLength(1) = 0.04;
% ax.YAxis.MinorTickValues = ax.YAxis.Limits(1):0.25:ax.YAxis.Limits(2);
% ax.XAxis.MinorTickValues = ax.XAxis.Limits(1):0.5:ax.XAxis.Limits(2);
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontsize',14)  % end inset
saveas(gcf, [folder_address '/', 'X.png'])
% close()