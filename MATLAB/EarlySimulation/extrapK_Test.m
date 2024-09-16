% Extrapolation of K in the limit of infinite lattice (zero inverse lattice size).
% Ks and Ls must be previously computed from generalcell.m and available in workspace 

%% plotting inverse system sizes vs dyadic elements

% inverse system sizes
invLs = 1./Ls;
m = size(Ks, 2);

hold on;
a = 1;

for i = 1:m
    for j = 1:m
        subplot(m,m,a)
        scatter(invLs, Ks(:,i,j),10 ,'filled')
        ylabel(['K_',num2str(i),'_',num2str(j)])
        xlabel('inverse length')
        
        a = a + 1;
    end
end
set(gcf, 'PaperUnits', 'inches');
 x_width=20 ;y_width=20;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
 saveas(gcf,'Triangular K vs invL.png')
%% finding extrapolated K elements


% inverse system sizes
invLs = 1./Ls;

% extrapolation of each element using linear method

m = size(Ks, 2);
Kextrap = zeros(m, m);

for i = 1:m
    for j = 1:m
        Kextrap(i, j) = interp1(invLs, Ks(:,i,j), 0, 'spline',  'extrap');
    end
end

global K_per_cell;
global K_per_dipole;
K_per_cell = Kextrap;
K_per_dipole = K_per_cell/n;



%% extrap methods test
clc

invLs = 1./Ls;
m = size(Ks, 2);

KextrapInterp = zeros(m, m);
KextrapPoly = zeros(m, m);

for i = 1:m
    for j = 1:m
        p = polyfit(invLs, Ks(:,i,j)', 3);
        KextrapPoly(i, j) = polyval(p, 0);
        KextrapInterp(i, j) = interp1(invLs, Ks(:,i,j), 0, 'linear',  'extrap');
    end
end

a = 1;

for i = 1:m
    for j = i:m
        subplot(m,m,m*(i-1) + j)
        hold on;
        scatter(invLs, Ks(:,i,j),3 ,'filled')
        scatter(0, KextrapPoly(i,j),3, 'red','filled')
        scatter(0, KextrapInterp(i,j),3, 'green','filled')
        ylabel(['K_',num2str(i),'_',num2str(j)])
%         xlabel('inverse length')
        
        a = a + 1;
    end
end
set(gcf, 'PaperUnits', 'inches');
 x_width=30 ;y_width=30;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
 saveas(gcf,'test.png')

