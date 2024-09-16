% Dyadic for Simple Unit Cell, 2D lattice with 3D dipoles

addpath(pwd);
addpath([pwd,'\lattice2d']);
clear all;
close all;
% clc;

tic();      %  start timer


% Constants 
% R1 and R2 are lattice vectors

% square lattice
% R1 = [1; 0];
% R2 = [0; 1];

% rectangular (1x2) lattice
% R1 = [1; 0]/(sqrt(2));
% R2 = [0; 2]/(sqrt(2));

% triangular lattice
R1 = [1; 0; 0];
R2 = [1/2; sqrt(3)/2; 0];

R = [R1 R2];
r = (R1 + R2)/2;

% Normalize vectors so that every lattice have density of one dipole per
% unit volume
% A = norm(cross(R1, R2));
% factor = nthroot(1/A, 2);
% 
% R = R * factor;
% r = r * factor;

% we want to find interaction of a dipole at (0, 0) with other dipoles 
% inside Radius L
L = 100;      

% total interaction per dipole
J_per_cell = JPerCellSimple(R, L)



% normalized dyadic
% J_normalized = J_per_cell/norm(J_per_cell)


toc()	 % stop timer
plotDiskLattice(R, r, 3)


% finding J for differnt lattice sizes to extrapolate J at infinite size
% Ls = [10 20 30 40 50 60 70 80 90 100 200 300 400 500 1000 ... 
%       ];
% Js = zeros(length(Ls), 3, 3);
% 
% for i = 1:length(Ls)
%     Js(i,:,:) = JPerCellCircleBnd(R1, R2, Ls(i));
% end










