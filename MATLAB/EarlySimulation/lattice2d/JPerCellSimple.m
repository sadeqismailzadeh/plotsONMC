function Jout = JPerCellSimple(R, L)
% compute total J dyadic per cell with a simple method
% total interaction of a dipole at (0, 0) with other dipoles in the lattice
% divided by 2. because an interaction is between 2 dipoles and each dipole
% contributes to half of this energy in total lattice Hamiltonian.
% 
% return:  J dyadic per dipole
% 
% R           3 * 2 matrix of lattice vectors   
% L           half size of lattice

R1 = R(:,1);
R2 = R(:,2);

Jout = zeros(3,3);

L2 = L*L;
for i = -L:L
    for j = -L:L
        if ((i == 0) && (j == 0))
            continue;
        end
        rij = i*R1 + j*R2;
        Jout = Jout + J(rij);
    end
end

Jout = Jout / 2;

end