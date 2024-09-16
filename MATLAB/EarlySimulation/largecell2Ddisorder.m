% Dyadic for Large Unit Cell, 2D lattice with 3D dipoles

% addpath(pwd);

% clear;
% close all;
% clc;

% tic();      %  start timer

% Constants

% R1 and R2 are lattice vectors
% r is 3*n matrix of Relative position of particles inside a cell, n is
% number of dipoles in a cell

  
% honeycomb lattice
% R1 = [1; 0; 0];
% R2 = [1/2; sqrt(3)/2; 0];
% r(:,1) = 1/3 * (R1 + R2);
% r(:,2) = 2/3 * (R1 + R2);



% rectangular (1x2) lattice
% R1 = [1; 0; 0];
% R2 = [0; 2; 0];
% r(:,1) = R1/2 +   R2/4;
% r(:,2) = R1/2 + 3*R2/4;


% % triangular lattice
R1 = [1; 0; 0];
R2 = [1/2; sqrt(3)/2; 0];
% r(:,1) = R1/2 +   R2/4;
% r(:,2) = R1/2 + 3*R2/4;
        

% square lattice
% R1 = [1; 0; 0];
% R2 = [0; 1; 0];
% r(:,1) = [0.25; 0.25; 0];
% r(:,2) = [0.25; 0.75; 0];
% r(:,3) = [0.75; 0.25; 0];
% r(:,4) = [0.75; 0.75; 0];


% large bravis lattice for arbitrary R1 and R2

r1 = (R1 + R2) / 2;

l = 4;
a = 1;

for i = 0:l-1
    for j = 0:l-1
        r(:,a) = (r1 + i*R1 + j*R2) / l;
        a = a + 1;
    end
end




% large honeycomb
% R1 = [1; 0; 0];
% R2 = [1/2; sqrt(3)/2; 0];
% r1 = 1/3 * (R1 + R2);
% r2 = 2/3 * (R1 + R2);
% 
% l = 4;
% a = 1;
% for i = 0:l-1
%     for j = 0:l-1
%         r(: ,a) = (r1 + i*R1 + j*R2) / l;
%         r(: ,a + 1) = (r2 + i*R1 + j*R2) / l;
%         a = a + 2;
%     end
% end



R = [R1 R2];

n = size(r,2);  % number of dipoles inside a cell
m = 3*n;        % auxiliary parameter

% Normalize vectors so that every lattice have density of one dipole per
% unit volume
A = norm(cross(R1, R2));
factor = sqrt(n/A);

R = R * factor;
r = r * factor;

% random positional disorder
% descriptions are in my notes
sigma_r = 0.1;
gaussMatrix = random('normal',0 ,sigma_r/sqrt(3) ,size(r))
% r = r + ((rand(3,size(r,2)) -0.5)*2) * 0.05;
r = r + gaussMatrix;

%% dyadic for fixed latttice  size

% we want to find interaction of a cell at (0, 0) with other cells 
% inside Radius L
% L = 10;          


% K_per_cell = KPerCellSimple(R, r ,L)
% K_per_dipole = K_per_cell/n

% Computing normalized J dyadic if all spins were identical
% J_sum = zeros(3,3);     
% for i = 1:n
%     for j = i:n
%         J_sum = J_sum + K_per_cell(3*i-2:3*i, 3*j-2:3*j);
%     end
% end
% J_per_dipole = J_sum/n
% J_normalized = J_per_dipole/norm(J_per_dipole)


% Dyadic K is normalized by dividing by its norm
% K_normalized = K_per_cell/norm(K_per_cell);
% K_normalized = round(K_normalized, 5)


%% extrapolate K and J

global K_per_cell;
global K_per_dipole;

Ls = [100  200 500 ... 
      ];

extrapKJ
% run('MCsimulation.m')
%% A visualization of lattice
plotLattice(R, r, 1)

% % 
% toc()	 % stop timer




