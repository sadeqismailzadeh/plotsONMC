function KpcUT = KPerCell(R, r, L)
% return: K dyadic per cell
% 
% interaction of dipoles inside a cell plus half of total interaction of 
% a cell at (0, 0) with other cells inside Radius L (cells inside a sphere 
% with radius L and center at (0,0)).
% we divide total interaction with other cells by 2. because an interaction is between 2 dipoles 
% and each dipole contributes to half of this energy in total lattice Hamiltonian.
%
% R        3 * 3 matrix of lattice vectors   
% r        3 * n  array of relettive particle positions inside a cell
% L        radius of sphere that bounds lattice

R1 = R(:,1);
R2 = R(:,2);
R3 = R(:,3);

n = size(r,2);  %number of dipoles inside a cell
m = 3*n;

% interaction of particles inside a cell
K_internal = zeros(m,m);    
for i = 1:n-1
    for j = i+1:n
        K_internal(3*i-2:3*i, 3*j-2:3*j) = J(r(:,i) - r(:,j));
    end
end


% total interaction of a cell at (0, 0) with other cells inside Radius L
K_tot = zeros(m,m); 

% this method works for any arbitrary combination of unit vectors
% but to understand what each part of code does, we assume:
% R1 = [1; 0; 0]
% R2 = [0; 1; 0]
% R3 = [0; 0; 1]

% accounts cells on x = 0 plane
i = 0;
j = 0;
[kmin, kmax] = kExtremes(R, L, i, j);
% kExtremes finds minimum and maximum value of k at fixed i and j such that 
% |i*R1 + j*R2 + k*R3| < L

while (isreal(kmin))
    for k = kmin:kmax
        if (i==0 && j==0 && k==0)
            continue;
        end
        rij = i*R1 + j*R2 + k*R3;
        K_tot = K_tot + K(rij, r, n);
    end
    
    j = j + 1;
    [kmin, kmax] = kExtremes(R, L, i, j);
end

i = 0;
j = -1;
[kmin, kmax] = kExtremes(R, L, i, j);
while (isreal(kmin))
    for k = kmin:kmax
        rij = i*R1 + j*R2 + k*R3;
        K_tot = K_tot + K(rij, r, n);
    end
    j = j - 1;
    [kmin, kmax] = kExtremes(R, L, i, j);
end



% accounts cells in first and eighth octants
% fisrt : (x, y, z) = (+, +, +)
% eighth: (x, y, z) = (+, +, -)
i = 1;
j = 0;
[kmin, kmax] = kExtremes(R, L, i, j);

while (isreal(kmin))
    while (isreal(kmin))
        for k = kmin:kmax
            rij = i*R1 + j*R2 + k*R3;
            K_tot = K_tot + K(rij, r, n);
        end
        j = j + 1;
        [kmin, kmax] = kExtremes(R, L, i, j);
    end
    i = i + 1;
    j = 0;
    [kmin, kmax] = kExtremes(R, L, i, j);
end



% accounts cells in fourth and fifth octants
% fourth: (x, y, z) = (+, -, +)
% fifth : (x, y, z) = (+, -, -)
i = 1;
j = -1;
[kmin, kmax] = kExtremes(R, L, i, j);

while (isreal(kmin))
    while (isreal(kmin))
        for k = kmin:kmax
            rij = i*R1 + j*R2 + k*R3;
            K_tot = K_tot + K(rij, r, n);
        end
        j = j - 1;
        [kmin, kmax] = kExtremes(R, L, i, j);
    end
    i = i + 1;
    j = -1;
    [kmin, kmax] = kExtremes(R, L, i, j);
end



% accounts cells in second and seventh octants
% second : (x, y, z) = (-, +, +)
% seventh: (x, y, z) = (-, +, -)
i = -1;
j = 0;
[kmin, kmax] = kExtremes(R, L, i, j);

while (isreal(kmin))
    while (isreal(kmin))
        for k = kmin:kmax
            rij = i*R1 + j*R2 + k*R3;
            K_tot = K_tot + K(rij, r, n);
        end
        j = j + 1;
        [kmin, kmax] = kExtremes(R, L, i, j);
    end
    i = i - 1;
    j = 0;
    [kmin, kmax] = kExtremes(R, L, i, j);
end



% accounts cells in third and sixth octants
% third: (x, y, z) = (-, -, +)
% sixth: (x, y, z) = (-, -, -)
i = -1;
j = -1;
[kmin, kmax] = kExtremes(R, L, i, j);

while (isreal(kmin))
    while (isreal(kmin))
        for k = kmin:kmax
            rij = i*R1 + j*R2 + k*R3;
            K_tot = K_tot + K(rij, r, n);
        end
        j = j - 1;
        [kmin, kmax] = kExtremes(R, L, i, j);
    end
    i = i - 1;
    j = -1;
    [kmin, kmax] = kExtremes(R, L, i, j);
end


% Total interaction per cell
Kpc = K_tot/2 + K_internal;


% With no effect on the Hamiltonian, Dyadic K can be made "Upper Traingular" 
% by Adding lower part to its upper counterpart
KpcUT = zeros(m,m);
for i = 1:n
    KpcUT(3*i-2:3*i, 3*i-2:3*i) = Kpc(3*i-2:3*i, 3*i-2:3*i);
end
for i = 1:n-1
    for j = i+1:n
        KpcUT(3*i-2:3*i, 3*j-2:3*j) = Kpc(3*i-2:3*i, 3*j-2:3*j) + ...
                                      Kpc(3*j-2:3*j, 3*i-2:3*i);
    end
end
end