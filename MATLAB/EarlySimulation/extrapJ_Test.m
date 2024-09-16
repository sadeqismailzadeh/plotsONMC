% Extrapolation of J in the limit of infinite lattice (zero inverse lattice size).
% Js and Ls must be previously computed from simple_cell.m and available in workspace 


%% plotting inverse system sizes vs dyadic elements

% inverse system sizes
invLs = 1./Ls;
hold on;
a = 1;
for i = 1:3
    for j = 1:3
        subplot(3,3,a)
        scatter(invLs, Js(:,i,j),15 ,'filled')
        ylabel(['J_',num2str(i),'_',num2str(j)])
        xlabel('inverse length')
        
        a = a + 1;
    end
end



%% finding extrapolated J elements


% inverse system sizes
invLs = 1./Ls;

% extrapolation of each element using linear method
Jextrap = zeros(3,3)
for i = 1:3
    for j = 1:3
        Jextrap(i, j) = interp1(invLs, Js(:,i,j), 0, 'linear',  'extrap');
    end
end


Jextrap








