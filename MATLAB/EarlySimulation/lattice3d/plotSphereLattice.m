function plotSphereLattice(R, r ,L)
% plot a visualization of lattice with spherical boundary
% 
% R      3 * 3  matrix of lattice vectors
% r      3 * n  matrix of relettive particle positions inside a cell
% l      radius of sphere

R1 = R(:,1);
R2 = R(:,2);
R3 = R(:,3);

n = size(r,2);  %number of dipoles inside a cell
a = 1;


% this method works for any arbitrary combination of unit vectors
% but to understand what each part of code does, we assume:
% R1 = [1; 0; 0]
% R2 = [0; 1; 0]
% R3 = [0; 0; 1]


% accounts cells on x = 0 plane
i = 0;
j = 0;
[kmin, kmax] = kExtremes(R, L, i, j);
% kExtremes finds minimum and maximum value of k at fixed i and j such that
% |i*R1 + j*R2 + k*R3| < L

while (isreal(kmin))
    for k = kmin:kmax
%         if (i==0 && j==0 && k==0)
%             continue;
%         end
        for t = 1: n
            x(a) = i*R1(1) + j*R2(1) + k*R3(1)+ r(1,t);
            y(a) = i*R1(2) + j*R2(2) + k*R3(2)+ r(2,t);
            z(a) = i*R1(3) + j*R2(3) + k*R3(3)+ r(3,t);
            
            a = a + 1;
        end
    end
    
    j = j + 1;
    [kmin, kmax] = kExtremes(R, L, i, j);
end


i = 0;
j = -1;
[kmin, kmax] = kExtremes(R, L, i, j);
while (isreal(kmin))
    for k = kmin:kmax
        for t = 1: n
            x(a) = i*R1(1) + j*R2(1) + k*R3(1)+ r(1,t);
            y(a) = i*R1(2) + j*R2(2) + k*R3(2)+ r(2,t);
            z(a) = i*R1(3) + j*R2(3) + k*R3(3)+ r(3,t);
            
            a = a + 1;
        end
    end
    
    j = j - 1;
    [kmin, kmax] = kExtremes(R, L, i, j);
end



% accounts cells in first and eighth octants
% fisrt : (x, y, z) = (+, +, +)
% eighth: (x, y, z) = (+, +, -)
i = 1;
j = 0;
[kmin, kmax] = kExtremes(R, L, i, j);

while (isreal(kmin))
    while (isreal(kmin))
        for k = kmin:kmax
            for t = 1: n
                x(a) = i*R1(1) + j*R2(1) + k*R3(1)+ r(1,t);
                y(a) = i*R1(2) + j*R2(2) + k*R3(2)+ r(2,t);
                z(a) = i*R1(3) + j*R2(3) + k*R3(3)+ r(3,t);
                
                a = a + 1;
            end
        end
        
        j = j + 1;
        [kmin, kmax] = kExtremes(R, L, i, j);
    end
    
    i = i + 1;
    j = 0;
    [kmin, kmax] = kExtremes(R, L, i, j);
end



% accounts cells in fourth and fifth octants
% fourth: (x, y, z) = (+, -, +)
% fifth : (x, y, z) = (+, -, -)
i = 1;
j = -1;
[kmin, kmax] = kExtremes(R, L, i, j);

while (isreal(kmin))
    while (isreal(kmin))
        for k = kmin:kmax
            for t = 1: n
                x(a) = i*R1(1) + j*R2(1) + k*R3(1)+ r(1,t);
                y(a) = i*R1(2) + j*R2(2) + k*R3(2)+ r(2,t);
                z(a) = i*R1(3) + j*R2(3) + k*R3(3)+ r(3,t);
                
                a = a + 1;
            end
        end
        
        j = j - 1;
        [kmin, kmax] = kExtremes(R, L, i, j);
    end
    
    i = i + 1;
    j = -1;
    [kmin, kmax] = kExtremes(R, L, i, j);
end


% accounts cells in second and seventh octants
% second : (x, y, z) = (-, +, +)
% seventh: (x, y, z) = (-, +, -)
i = -1;
j = 0;
[kmin, kmax] = kExtremes(R, L, i, j);

while (isreal(kmin))
    while (isreal(kmin))
        for k = kmin:kmax
            for t = 1: n
                x(a) = i*R1(1) + j*R2(1) + k*R3(1)+ r(1,t);
                y(a) = i*R1(2) + j*R2(2) + k*R3(2)+ r(2,t);
                z(a) = i*R1(3) + j*R2(3) + k*R3(3)+ r(3,t);
                
                a = a + 1;
            end
        end
        
        j = j + 1;
        [kmin, kmax] = kExtremes(R, L, i, j);
    end
    
    i = i - 1;
    j = 0;
    [kmin, kmax] = kExtremes(R, L, i, j);
end



% accounts cells in third and sixth octants
% third: (x, y, z) = (-, -, +)
% sixth: (x, y, z) = (-, -, -)
i = -1;
j = -1;
[kmin, kmax] = kExtremes(R, L, i, j);

while (isreal(kmin))
    while (isreal(kmin))
        for k = kmin:kmax
            for t = 1: n
                x(a) = i*R1(1) + j*R2(1) + k*R3(1)+ r(1,t);
                y(a) = i*R1(2) + j*R2(2) + k*R3(2)+ r(2,t);
                z(a) = i*R1(3) + j*R2(3) + k*R3(3)+ r(3,t);
                
                a = a + 1;
            end
        end
        
        j = j - 1;
        [kmin, kmax] = kExtremes(R, L, i, j);
    end
    
    i = i - 1;
    j = -1;
    [kmin, kmax] = kExtremes(R, L, i, j);
end


scatter3(x, y, z,'filled');
axis('equal');


end