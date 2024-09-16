function Jout = JPerCellSphereBnd(R, L)
% total interaction of a dipole at (0, 0) with other dipoles inside Radius L
% (dipoles inside a Sphere with radius L and center at (0,0))
% this value is divided by 2. because an interaction is between 2 dipoles 
% and each dipole contributes to half of this energy in total lattice Hamiltonian.
%
% return:  J dyadic per dipole
%
% R      3 * 3  matrix of lattice vectors 
% L      half size of lattice

R1 = R(:,1);
R2 = R(:,2);
R3 = R(:,3);

Jout = zeros(3,3);


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
        Jout = Jout + J(rij);
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
        Jout = Jout + J(rij);
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
            Jout = Jout + J(rij);
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
            Jout = Jout + J(rij);
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
            Jout = Jout + J(rij);
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
            Jout = Jout + J(rij);
        end
        j = j - 1;
        [kmin, kmax] = kExtremes(R, L, i, j);
    end
    i = i - 1;
    j = -1;
    [kmin, kmax] = kExtremes(R, L, i, j);
end


% Jpercell must be divided by 2
Jout = Jout / 2;

end
