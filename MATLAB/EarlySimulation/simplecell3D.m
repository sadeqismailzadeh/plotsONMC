% Dyadic for Simple Unit Cell, 3D lattice

addpath(pwd);
addpath([pwd,'\lattice3d']);
clear all;
close all;
clc;

tic();      %  start timer


% Constants 
% R1 and R2 are lattice vectors

% simple cubic
% R1 = [1; 0; 0];
% R2 = [0; 1; 0];
% R3 = [0; 0; 1];

% bcc
% R1 = [1; 0; 0];
% R2 = [0; 1; 0];
% R3 = [0.5; 0.5; 0.5];

% fcc
% R1 = [1; 1; 0];
% R2 = [1; 0; 1];
% R3 = [0; 1; 1];

R3 = [0; 1; 0];
R1 = [1; 0; 0];
R2 = [1/2; sqrt(3)/2; 0];


R = [R1 R2 R3];

% Normalize vectors so that every lattice have density of one dipole per
% unit volume
% V = abs(det(R));
% factor = nthroot(1/V, 3);
% R = R * factor;

% we want to find interaction of a dipole at (0, 0) with other dipoles 
% inside Radius L
L = 60;      
% total interaction per dipole
J_per_cell = JPerCellSimple(R, L)
% J_per_cell = JPerCellSimple(R, L)


% [kmin, kmax] = kExtremes (R, L, 60, 80)

% normalized dyadic
% J_normalized = J_per_cell/norm(J_per_cell)



% finding J for differnt lattice sizes to extrapolate J at infinite size
% Ls = [10 20 30 40 50 60 70 80 90 100 200 300 400 500 1000 ... 
%       ];
% Js = zeros(length(Ls), 3, 3);
% 
% for i = 1:length(Ls)
%     Js(i,:,:) = JPerCellCircleBnd(R1, R2, Ls(i));
% end

% plotLattice(R, 3);

toc()	 % stop timer





