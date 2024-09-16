function h = Hamiltonian(A)
%HAMILTONIAN Hamiltonian of Lattice per dipole
% assuming magnitude of each dipole is equal to 1
% A     spherical angles of dipoles, A(2i-1) -> theta(i), A(2i) -> phi(i)

n = round(length(A)/2);
h = 0;

global K_per_dipole;

for i = 1:n
    for j = i:n 
        h = h - ... 
            [sin(A(2*i-1))*cos(A(2*i));
             sin(A(2*i-1))*sin(A(2*i));
             cos(A(2*i-1))]' ...     
             * K_per_dipole(3*i-2:3*i, 3*j-2:3*j) * ...    
             [sin(A(2*j-1))*cos(A(2*j));
              sin(A(2*j-1))*sin(A(2*j));
              cos(A(2*j-1))];
    end
end
            
end

