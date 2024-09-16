function dipoles = run1step(dipoles, KpcSym, T)
% run one montecarlo step and changes state of dipoles
% input:
%   dipoles       3 * n stores dipole moments' vectors.
%                 i.e M(5,1) represents magnitudeof 5th dipole along x direction.
%   KpcSym        symmetric version of Kpercell dyadic
%   T             Temperature of configuration
%
% output:
%   #changed version of dipoles

n = length(dipoles) / 3;
for c = 1 : n
    i = randi(n);
    
%     li = norm(Moments(:,i));
    li = 1;
%     mi_new = Moments(:,i) + 0.1*(rand(3,1) - 0.5);
%     mi_new = [rand(2,1) - 0.5; 0];
    mi_new = rand(3,1) -0.5;
    mi_new = mi_new * li / norm(mi_new);
    
    dE = dErot(dipoles, i, mi_new, KpcSym, n);
    
    if dE < 0
        dipoles(3*i-2: 3*i) = mi_new;
    elseif (dE > 0) && (rand < exp(-dE / T))
        dipoles(3*i-2: 3*i) = mi_new;
    end
end

