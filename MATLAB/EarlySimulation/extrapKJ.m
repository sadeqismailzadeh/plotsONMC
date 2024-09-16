% Computation of dyadic in different lattice sizes and extrapolate at infinite latttice size
% depends on largecell2(3)D.m

global K_per_cell;
global K_per_dipole;

Js = zeros(length(Ls), 3, 3);
Ks = zeros(length(Ls), m, m);

for k = 1:length(Ls)
    Kpc = KPerCell(R, r, Ls(k));
    Ks(k, :, :) = Kpc;

    J_sum = zeros(3,3);     
    for i = 1:n
        for j = i:n
            J_sum = J_sum + Kpc(3*i-2:3*i, 3*j-2:3*j);
        end
    end
    Js(k, :, :) = J_sum/n;
    L = Ls(k)
end

invLs = 1./Ls;
Kextrap = zeros(m, m);
for i = 1:m
    for j = 1:m
        Kextrap(i, j) = interp1(invLs, Ks(:,i,j), 0, 'spline',  'extrap');
    end
end
K_per_cell = Kextrap;
K_per_dipole = K_per_cell/n;

% Computing normalized J dyadic if all spins were identical
Jextrap = zeros(3,3);
for i = 1:3
    for j = 1:3
        Jextrap(i, j) = interp1(invLs, Js(:,i,j), 0, 'spline',  'extrap');
    end
end
Jextrap

% with no change in Hamiltonian we can make K dyadic a symmetric matrix symmetric by dividing
% non-diagonal elements between upper and lower part (Note: dyadic is already upper triangular).
% a symmetric matrix shows better results in eigenvalue problem. because it is guaranteed to have
% real eigenvalues and its eigenvectors form a complete set of mutually orthonormal basis
KpcSym = K_per_cell;
for i = 1:n
    for j = i:n
        KpcSym(3*j-2:3*j, 3*i-2:3*i) = KpcSym(3*j-2:3*j, 3*i-2:3*i) + KpcSym(3*i-2:3*i, 3*j-2:3*j);
    end
end
KpcSym = KpcSym / 2;
KpcUT = K_per_cell;

