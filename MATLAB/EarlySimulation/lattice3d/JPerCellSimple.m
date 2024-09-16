function Jout = JPerCellSimple(R, L)
% total interaction of a dipole at (0, 0) with other dipoles in the lattice
% divided by 2. because an interaction is between 2 dipoles and each dipole
% contributes to half of this energy in total lattice Hamiltonian.

% return:  J dyadic per dipole

% R           lattice vectors
% L           half size of lattice

R1 = R(:,1);
R2 = R(:,2);
R3 = R(:,3);

Jout = zeros(3,3);

L2 = L*L;
for i = -L:L
    for j = -L:L
        for k = -L:L
            if ((i == 0) && (j == 0) && (k == 0))
                continue;
            end
            rij = i*R1 + j*R2 + k*R3;
            Jout = Jout + J(rij);
        end
    end
end

Jout = Jout / 2;

end
