function plotDiskLattice(R, r, L)
% plot a visualization of lattice with circular boundary
% 
% R           3 * 2 matrix oflattice vectors
% r           3 * n  array of relettive particle positions inside a cell
% L           radius of circle that bounds lattice

R1 = R(:,1);
R2 = R(:,2);

n = size(r,2);  %number of dipoles inside a cell
a = 1;

% this method works for any arbitrary combination of unit vectors
% but to understand what each part of code does, we assume:
% R1 = [1; 0; 0]
% R2 = [0; 1; 0]

% accounts cells on x = 0 line
i = 0;
% jExtremes finds minimum and maximum value of j at fixed i such that 
% |i*R1 + j*R2| < L
[jmin, jmax] = jExtremes(R1, R2, L, i);

for j = jmin:jmax
    if (i==0 && j==0)
        continue;
    end 
    for k = 1: n
        x(a) = i*R1(1)+j*R2(1)+r(1,k);
        y(a) = i*R1(2)+j*R2(2)+r(2,k);
    
        a = a + 1;
    end
end

% accounts cells in first and fourth quadrants
% fisrt : (x, y) = (+, +)
% fourth: (x, y) = (+, -)
i = 1;
[jmin, jmax] = jExtremes(R1, R2, L, i);

while (isreal(jmin))
    for j = jmin:jmax
        for k = 1: n
            x(a) = i*R1(1)+j*R2(1)+r(1,k);
            y(a) = i*R1(2)+j*R2(2)+r(2,k);
                
            a = a + 1;
        end
    end
    
    i = i + 1;
    [jmin, jmax] = jExtremes(R1, R2, L, i);
end

% accounts cells in second and third quadrants
% second: (x, y) = (-, +)
% third : (x, y) = (-, -)
i = -1;
[jmin, jmax] = jExtremes(R1, R2, L, i);

while (isreal(jmin))
    for j = jmin:jmax
        for k = 1: n
            x(a) = i*R1(1)+j*R2(1)+r(1,k);
            y(a) = i*R1(2)+j*R2(2)+r(2,k);
                
            a = a + 1;
        end
    end
    
    i = i - 1;
    [jmin, jmax] = jExtremes(R1, R2, L, i);
end

scatter(x, y,'filled');
axis('equal');


end