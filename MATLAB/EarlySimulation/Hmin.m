% Minimization of Hamiltonian

latticeName = 'square';
clc

tic();

ni = 2*n;


Ai0 = rand(1,ni)*2*pi;

[Ai, E_per_dip] = fminunc(@Hamiltonian,Ai0);

A = zeros(n, 2);

for i = 1:n
    A(i, 1) = Ai(2*i-1);
    A(i, 2) = Ai(2*i);
end

A*180/pi

for i = 1:n
    A(i,1) = mod(A(i,1), 2*pi);
    A(i,2) = mod(A(i,2), 2*pi);
    
    if  (A(i,1) > pi) && (A(i,1) < 2*pi)
        A(i,1) = 2*pi - A(i,1);
        A(i,2) = mod(A(i,2)+pi, 2*pi);
    end
end

theta___phi = A*180/pi
E_per_dip

l = ones(n,1) * 0.15;

toc()
%% save
addpath(pwd);
addpath([pwd,'\Lattice3D']);
plotLatticeDipole(R, r, 2, A, l);
% set(gcf, 'PaperUnits', 'inches');
%  x_width=10 ;y_width=10;
%  set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf,[latticeName,', n= ', num2str(n),', E=', num2str(E_per_dip),' .png'])







