function dE = dErot(dipoles, i, mi_new, KpcSym, n)
% energy difference of ith dipole by converting to mi_new
% input:
% dipoles       3 * n stores dipole moments' vectors.
%               i.e M(5,1) represents magnitudeof 5th dipole along x direction.
% i             index of changing dipole
% mi_new        new state of ith dipole
% KpcSym        symmetric version of Kpercell dyadic
%n              number of dipoles in one cell

dE = -(mi_new - dipoles(3*i-2: 3*i))' * KpcSym(3*i-2:3*i, 3*i-2:3*i) * (mi_new + dipoles(3*i-2: 3*i));
   
dM = mi_new - dipoles(3*i-2: 3*i);

for j = 1:n
    if j ~= i
        dE = dE - 2*dM'*KpcSym(3*i-2:3*i, 3*j-2:3*j)*dipoles(3*j-2: 3*j);
    end
end

end

