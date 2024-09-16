function dipoles = stabilize(dipoles, KpcSym, T, Nsteps)
% stabilize state of dipoles after changing temperature
% input:
%   dipoles       3 * n stores dipole moments' vectors.
%                 i.e M(5,1) represents magnitudeof 5th dipole along x direction.
%   KpcSym        symmetric version of Kpercell dyadic
%   Nsteps        number of stabilization steps
%
% output:
%   #changed version of dipoles

for c = 1: Nsteps
    dipoles = run1step(dipoles, KpcSym, T);
end

end