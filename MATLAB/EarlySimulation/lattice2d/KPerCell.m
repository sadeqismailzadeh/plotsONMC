function KpcUT = KPerCell(R, r, L)
% return: K dyadic per cell
% 
% interaction of dipoles inside a cell plus half of total interaction of 
% a cell at (0, 0) with other cells inside Radius L (cells inside a disk 
% with radius L and center at (0,0)).
% we divide total interaction with other cells by 2. because an interaction is between 2 dipoles 
% and each dipole contributes to half of this energy in total lattice Hamiltonian.
%
% R        3 * 2 matrix of lattice vectors   
% r        3 * n  array of relettive particle positions inside a cell
% L        radius of circle that bounds lattice

R1 = R(:,1);
R2 = R(:,2);

n = size(r,2);  % number of dipoles inside a cell
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

% accounts cells on x = 0 line
i = 0;
% jExtremes finds minimum and maximum value of j at fixed i such that 
% |i*R1 + j*R2| < L
[jmin, jmax] = jExtremes(R1, R2, L, i);

for j = jmin:jmax
    if (i==0 && j==0)
        continue;
    end
    rij = i*R1 + j*R2;
    K_tot = K_tot + K(rij, r, n);
end

% accounts cells in first and fourth quadrants
% fisrt : (x, y) = (+, +)
% fourth: (x, y) = (+, -)
i = 1;
[jmin, jmax] = jExtremes(R1, R2, L, i);

while (isreal(jmin))
    for j = jmin:jmax
        rij = i*R1 + j*R2;
        K_tot = K_tot + K(rij, r, n);
    end
    
    i = i + 1;
    [jmin, jmax] = jExtremes(R1, R2, L, i);
end

% accounts cells in second and third quadrants
% second: (x, y) = (-, +)
% third : (x, y) = (-, -)
i = -1;
[jmin, jmax] = jExtremes(R1, R2, L, i);

while (isreal(jmin))
    for j = jmin:jmax
        rij = i*R1 + j*R2;
        K_tot = K_tot + K(rij, r, n);
    end
    
    i = i - 1;
    [jmin, jmax] = jExtremes(R1, R2, L, i);
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
