function  [jmin, jmax] = jExtremes(R1, R2, L, i)
% finds minimum and maximum value of j at fixed i such that |i*R1 + j*R2| < L

R1R2 = R1'*R2;  
R12 = R1'*R1;
R22 = R2'*R2;

sqrdelta = sqrt(4*i*i*(R1R2)^2 - 4*(R22* (i*i*R12 - L^2)));

j1 = (-2*i*R1R2 + sqrdelta)/(2*R22);
j2 = (-2*i*R1R2 - sqrdelta)/(2*R22);

if j1 > j2
    jmax = j1;
    jmin = j2;
else
    jmax = j2;
    jmin = j1;
end

jmax = floor(jmax);
jmin = ceil(jmin);

end
