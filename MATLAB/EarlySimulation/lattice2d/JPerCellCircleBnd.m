function Jout = JPerCellCircleBnd(R, L)
% compute total J dyadic per cell with circular lattice boundary
% total interaction of a dipole at (0, 0) with other dipoles inside Radius L
% (dipoles inside a disk with radius L and center at (0,0))
% this value is divided by 2. because an interaction is between 2 dipoles 
% and each dipole contributes to half of this energy in total lattice Hamiltonian.
% 
% return:  J dyadic per dipole
% 
% R           3 * 2 matrix of lattice vectors    
% L           radius of circle that bounds lattice

R1 = R(:,1);
R2 = R(:,2);

% this method works for any arbitrary combination of unit vectors
% but to understand what each part of code does, we assume:
% R1 = [1; 0; 0]
% R2 = [0; 1; 0]

% accounts cells on x = 0 line
Jout = zeros(3,3);

i = 0;
% jExtremes finds minimum and maximum value of j at fixed i such that 
% |i*R1 + j*R2| < L
[jmin, jmax] = jExtremes(R1, R2, L, i);

for j = jmin:jmax
    if (i==0 && j==0)
        continue;
    end
    rij = i*R1 + j*R2;
    Jout = Jout + J(rij);
end

% accounts cells in first and fourth quadrants
% fisrt : (x, y) = (+, +)
% fourth: (x, y) = (+, -)
i = 1;
[jmin, jmax] = jExtremes(R1, R2, L, i);

while (isreal(jmin))
    for j = jmin:jmax
        rij = i*R1 + j*R2;
        Jout = Jout + J(rij);
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
        Jout = Jout + J(rij);
    end
    
    i = i - 1;
    [jmin, jmax] = jExtremes(R1, R2, L, i);
end


% Jpercell must be divided by 2
Jout = Jout / 2;

end