%% Minimization of Hamiltonian using largest eigenvectors of dyadic
clc
addpath(pwd);
addpath([pwd,'\lattice2d']);
tic()
% etot = eig(K_symmetric)
% K_symmetric1 = round(KpcSym, 3);
% e = eigs(KpcSym, 3, 'la')
[V ,~] = eigs(KpcSym, 6, 'la');
es = eig(KpcSym)/2

% V(:, 1) + 2*V(: ,2)

toc()

%% Check Hamiltonian
clc
tic()
H = @(M) - M'* KpcSym * M;
% latticeName = 'water 1.5';
% 
% a(1) = 1;
% a(2) = -1;
% a = a/norm(a);
% M1 = a(1)*V(:,1) + a(2)*V(:,2);
% V1 = V(:,1);
% V2 = V(:,2);
% H(V1) - H(V2)
% -(V1-V2)' * KpcSym * (V1 + V2)

% for Vindex = 1:size(V,2)
% M = V(:,Vindex);
M = V(:,1)

M = M * sqrt(n);

l = zeros(n ,1);
A = zeros(n, 2);

for i = 1:n
    x = M(3*i - 2);
    y = M(3*i - 1);
    z = M(3*i);
    
    [A(i,2), A(i,1), l(i)] = cart2sph(x, y, z);
    % second output is angle of elevation from xy plane
    
    A(i, 1) = pi/2 - A(i, 1); % convert to common spherical theta angle
    
end

M_per_dip = zeros(3,1);
Mx = 0;
My = 0;
Mz = 0;
for i = 1:n
    Mx = Mx + M(3*i - 2);
    My = My + M(3*i - 1);
    Mz = Mz + M(3*i);
end
    
[M_per_dip(3), M_per_dip(2), M_per_dip(1)] = cart2sph(Mx, My, Mz);

M_per_dip(2) = pi/2 - M_per_dip(2);
M_per_dip(1) = M_per_dip(1) / n;



r__theta__phi = [l, A*180/pi]

E_per_dip = round(H(M), 3, 'significant')
M_per_dip = round(M_per_dip, 3)

l = l*0.15;



clf
plotLatticeDipole(R, r, 1, A, l);
% view(40,20);
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% filename = [latticeName,', n= ', num2str(n),', config(',num2str(Vindex), ...
%             '), E=', num2str(E_per_dip),', M= ', num2str(M_per_dip(1))];
%     
% saveas(gcf,[filename,' .png'])
% saveas(gcf,[filename,' .fig'])
% end
% close all;
toc()
%% save
axis off
saveas(gcf, ['max_eigenvector.png'])
saveas(gcf, ['max_eigenvector.fig'])
close();

%% water
% clc
% addpath(pwd);
% addpath([pwd,'\lattice3d']);
% 
% H = @(M) - M'* KpcSym * M;
% 
% A = [0      0;
%      90    30;
%      90    150;
%      90    90;
%      0      0;
%      90     -90;
%      90     -30;
%      90     210;] * pi/180;
%  
%  l = ones(8,1);
%  M = zeros(3*n,1);
%  for i = 1:n
%     M(3*i - 2) = l(i) * sin(A(i,1)) * cos(A(i,2));
%     M(3*i - 1) = l(i) * sin(A(i,1)) * sin(A(i,2));
%     M(3*i    ) = l(i) * cos(A(i,1));
%  end
% 
%  H(M)
%  
%  l = l*0.15;
% plotLatticeDipole(R, r, 2, A, l);











