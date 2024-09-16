% Dyadic for Large Unit Cell, 3D lattice

addpath(pwd);
addpath([pwd,'\Lattice3D']);
clear all;
close all;
clc;

tic();      %  start timer

% Constants

% R1, R2, R3 are lattice vectors

% r is 3*n matrix of Relative position of particles inside a cell, n is
% number of dipoles in a cell

% honeycomb lattice
% R1 = [1; 0; 0];
% R2 = [1/2; sqrt(3)/2; 0];
% R3 = [0; 0; 2];
% 
% r(:,1) = 1/3 * (R1 + R2);
% r(:,2) = 2/3 * (R1 + R2);
% r(:,3) = 1/3 * (R1 + R2) + R3/2;
% r(:,4) = 2/3 * (R1 + R2) + R3/2;


% simple cubic
R1 = [1; 0; 0];
R2 = [0; 1; 0];
R3 = [0; 0; 1];

% r(:,1) =   R1/4 +   R2/4 +   R3/4;
% r(:,2) = 3*R1/4 +   R2/4 +   R3/4;  
% r(:,3) =   R1/4 + 3*R2/4 +   R3/4;
% r(:,4) =   R1/4 +   R2/4 + 3*R3/4;
% r(:,5) = 3*R1/4 + 3*R2/4 +   R3/4;
% r(:,6) = 3*R1/4 +   R2/4 + 3*R3/4;
% r(:,7) =   R1/4 + 3*R2/4 + 3*R3/4;
% r(:,8) = 3*R1/4 + 3*R2/4 + 3*R3/4;

% fcc 
% R1 = [1; 1; 0];
% R2 = [1; 0; 1];
% R3 = [0; 1; 1];


% bcc
% R1 = [1; 0; 0];
% R2 = [0; 1; 0];
% R3 = [0.5; 0.5; 0.5];


% large bravis lattice
% arbitrary R1 and R2

r1 = (R1 + R2) / 2;

l = 2;
a = 1;

for i = 0:l-1
    for j = 0:l-1
        for k =0:l-1
            r(: ,a) = (r1 + i*R1 + j*R2 + k*R3) / l;
            a = a + 1;
        end
    end
end
% r = r + rand(3,size(r,2)) / (4*l);

% water
% R1 = [1; 0; 0] / (3^0.5/4);
% R2 = [1/2; sqrt(3)/2; 0] / (3^0.5/4);
% R3 = [0; 0; 10];
% ri = 1/3 * (R1 + R2);
% rj = 2/3 * (R1 + R2);
% 
% b = norm(rj - ri)/2;
% r(: ,1) = ri + R3/4;
% r(: ,2) = ri + b*[cosd(30); sind(30); 0] + R3/2  ;
% r(: ,3) = ri + b*[cosd(150); sind(150); 0] + R3/2;
% r(: ,4) = ri + b*[cosd(-90); sind(-90); 0]+ R3/2 ;
% 
% r(: ,5) = rj + 3*R3/4;
% r(: ,6) = rj + b*[cosd(90); sind(90); 0]  + R3 ;
% r(: ,7) = rj + b*[cosd(-30); sind(-30); 0] + R3;
% r(: ,8) = rj + b*[cosd(210); sind(210); 0]+ R3 ;




% l = 1;
% a = 1;
% for i = 0:l-1
%     for j = 0:l-1
%         r(: ,a) =     (r1 + i*R1 + j*R2) + R3/4;
%         r(: ,a + 1) = (r2 + i*R1 + j*R2) + R3/2;
%         r(: ,a + 2) = (r3 + i*R1 + j*R2) + R3/2;
%         r(: ,a + 3) = (r4 + i*R1 + j*R2) + R3/2;
%         r(: ,a + 4) = (r5 + i*R1 + j*R2) + 3*R3/4;
%         r(: ,a + 5) = (r6 + i*R1 + j*R2) + R3;
%         r(: ,a + 6) = (r7 + i*R1 + j*R2) + R3;
%         r(: ,a + 7) = (r8 + i*R1 + j*R2) + R3;
%         
%         a = a + 8;
%     end
% end
% 
% R1 = R1 * l;
% R2 = R2 * l;


R = [R1 R2 R3];


% r(:,1) = R1/2 +   R2/4;
% r(:,2) = R1/2 + 3*R2/4;


n = size(r,2);  % number of dipoles inside a cell
m = 3*n;        % auxiliary parameter


% Normalize vectors so that every lattice have density of one dipole per
% unit volume
V = abs(det(R));
factor = nthroot(n/V, 3);

R = R * factor;
r = r * factor;

% A = norm(cross([R1], [R2]));
% factor = sqrt(n/A);
% 
% R = R * factor;
% r = r * factor;

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

Ls = [10 15 20 25 30 35 40   ... 
      ];

run('extrapKJ.m')
%% A visualization of lattice

plotCell(R ,r)
% plotSphereLattice(R, r, 2.5)
% plotLattice(R, r, 2)

toc()	 % stop timer

































