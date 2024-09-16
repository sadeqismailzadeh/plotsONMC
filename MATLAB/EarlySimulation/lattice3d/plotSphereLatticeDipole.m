function [ output_args ] = plotSphereLatticeDipole(R, r, L, A, l)
%PLOTSPHERELATTICEDIPOLE plots a visualization of spherical shaped lattice including dipoles
% R        3 * n matrix of lattice vectors
% r        3 * n  array of relettive particle positions inside a cell
% A        array of spherical Angles of dipoles,
%          theta of ith  dipole = A(1,i), phi of ith dipole = A(2,i)
% l        length of dipoles (all lenghts are assumed to be same)
% L        lattice size

R1 = R(:,1);
R2 = R(:,2);
R3 = R(:,3);

n = size(r,2);  %number of dipoles inside a cell


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
            x = i*R1(1) + j*R2(1) + k*R3(1)+ r(1,t);
            y = i*R1(2) + j*R2(2) + k*R3(2)+ r(2,t);
            z = i*R1(3) + j*R2(3) + k*R3(3)+ r(3,t);
            
            plotDipole([x, y, z], A(t, :), l(t));
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
            x = i*R1(1) + j*R2(1) + k*R3(1)+ r(1,t);
            y = i*R1(2) + j*R2(2) + k*R3(2)+ r(2,t);
            z = i*R1(3) + j*R2(3) + k*R3(3)+ r(3,t);
            
            plotDipole([x, y, z], A(t, :), l(t));
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
                x = i*R1(1) + j*R2(1) + k*R3(1)+ r(1,t);
                y = i*R1(2) + j*R2(2) + k*R3(2)+ r(2,t);
                z = i*R1(3) + j*R2(3) + k*R3(3)+ r(3,t);
                
                plotDipole([x, y, z], A(t, :), l(t));
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
                x = i*R1(1) + j*R2(1) + k*R3(1)+ r(1,t);
                y = i*R1(2) + j*R2(2) + k*R3(2)+ r(2,t);
                z = i*R1(3) + j*R2(3) + k*R3(3)+ r(3,t);
                
                plotDipole([x, y, z], A(t, :), l(t));
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
                x = i*R1(1) + j*R2(1) + k*R3(1)+ r(1,t);
                y = i*R1(2) + j*R2(2) + k*R3(2)+ r(2,t);
                z = i*R1(3) + j*R2(3) + k*R3(3)+ r(3,t);
                
                plotDipole([x, y, z], A(t, :), l(t));
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
                x = i*R1(1) + j*R2(1) + k*R3(1)+ r(1,t);
                y = i*R1(2) + j*R2(2) + k*R3(2)+ r(2,t);
                z = i*R1(3) + j*R2(3) + k*R3(3)+ r(3,t);
                
                plotDipole([x, y, z], A(t, :), l(t));
            end
        end
        
        j = j - 1;
        [kmin, kmax] = kExtremes(R, L, i, j);
    end
    
    i = i - 1;
    j = -1;
    [kmin, kmax] = kExtremes(R, L, i, j);
end

% xlim([-L-R1(1)+R2(1)+R3(1)  L+R1(1)+R2(1)+R3(1)]);
% ylim([-L-R1(2)+R2(2)+R3(2)  L+R1(2)+R2(2)+R3(2)]);
% zlim([-L-R1(3)+R2(3)+R3(3)  L+R1(3)+R2(3)+R3(3)]);
axis equal


end

