function [dipoles, meanE, C, orderparam, X] = getStatisticalData(dipoles, KpcSym, KpcUT, ...
                                                T, Nsteps, multiplier_to_decorrelate)
% compute statistical data by averaging over many monte carlo steps
% input:
%   dipoles       3 * n stores dipole moments' vectors.
%                 i.e M(5,1) represents magnitudeof 5th dipole along x direction.
%   KpcSym        symmetric version of Kpercell dyadic
%   KpcUt         upper triangular version of Kpercell dyadic
%   T             Temperature of configuration
%   Nsteps        number of averaging steps
%
% output:
%   dipoles       #changed version of dipoles
%   meanE         average energy
%   C             heat capacity
%   orderparam    order parameter
%   X             susceptibility 

n = length(KpcSym) / 3;
sumE = 0;
sumE2 = 0;
orderparam = 0;
orderparam2 = 0;
global maxeigvecs;
global numvecs;

for c = 1: Nsteps
    dipoles = run1step(dipoles, KpcSym, T);
    if mod(Nsteps, multiplier_to_decorrelate) == 0
        E_config = 0;
        orderparam_config = 0;


        E_config = - dipoles' * KpcUT * dipoles;

        for v = 1:numvecs
            orderparam_config = orderparam_config + (dot(dipoles, maxeigvecs(:, v))) ^2;
        end

        orderparam = orderparam + orderparam_config;
        orderparam2 = orderparam2 + orderparam_config^2;
        sumE = sumE + E_config;
        sumE2 = sumE2 + E_config^2;
    end
end

orderparam = orderparam /(Nsteps); %it shouldnt be Nsteps, consider decorrelate
orderparam2 = orderparam2 /(Nsteps);
meanE = sumE / (Nsteps);
meanE2 = sumE2 / (Nsteps);

% per dipole
X = (orderparam2 - orderparam^2) / (n * T);
C = (meanE2 - meanE ^2) / (n * T^2);
meanE = meanE/n;
orderparam = orderparam/n;
end

