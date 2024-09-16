function [cnum, dmin] = coordnum(R, r )
%COORDNUM Summary of this function goes here
%   Detailed explanation goes here


R1 = R(:,1);
R2 = R(:,2);
n = size(r,2);

rs = zeros(3, 9*n);
a = 1;

for i = -1:1
    for j = -1:1
        for t = 1:n
            if ( i==0 && j==0 && t==1)
                continue;
            end
            
            rs(:, a) = i*R1 + j*R2 + r(:, t);
            a = a + 1;
        end
    end
end


ds = zeros(a-1, 1);
for i = 1:a-1   
    ds(i) = norm(rs(:,i) - r(:,1));
end

tol = 0.0001;
cnum = length( find((ds > min(ds)-tol)  &  (ds < min(ds)+tol)) );
dmin = min(ds);

end

