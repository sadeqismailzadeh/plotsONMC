function  [kmin, kmax] = kExtremes(R, L, i, j)
% finds minimum and maximum value of k at fixed i and j such that 
% |i*R1 + j*R2 + k*R3| < L

a = R(:,3)'*R(:,3);
b = 2*i*R(:,1)'*R(:,3) + 2*j*R(:,2)'*R(:,3);
c = i*i*R(:,1)'*R(:,1) + j*j*R(:,2)'*R(:,2) + 2*i*j*R(:,1)'*R(:,2) - L^2;

sqrdelta = sqrt(b^2 - 4*a*c);

k1 = (-b + sqrdelta)/(2*a);
k2 = (-b - sqrdelta)/(2*a);

if k1 > k2
    kmax = k1;
    kmin = k2;
else
    kmax = k2;
    kmin = k1;
end

kmax = floor(kmax);
kmin = ceil(kmin);

end
