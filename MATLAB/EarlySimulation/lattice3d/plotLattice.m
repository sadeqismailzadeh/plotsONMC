function plotLattice(R, r, l)
% plot a visualization of lattice including cell boundaries
%
% R      3 * 3  matrix of lattice vectors
% r      3 * n  matrix of relettive particle positions inside a cell
% l      size of lattice

R1 = R(:,1);
R2 = R(:,2);
R3 = R(:,3);

N = l*l;
n = size(r,2);  %number of dipoles inside a cell

x = zeros(1, n*N);
y = zeros(1, n*N);
z = zeros(1, n*N);

% storing location of dipoles in x and y arrays
a = 1;
for i = 0:l-1
    for j = 0:l-1
        for k = 0:l-1
            
            for t = 1: n
                x(a) = i*R1(1) + j*R2(1) + k*R3(1)+ r(1,t);
                y(a) = i*R1(2) + j*R2(2) + k*R3(2)+ r(2,t);
                z(a) = i*R1(3) + j*R2(3) + k*R3(3)+ r(3,t);
                
                a = a + 1;
            end
            
        end
    end
end

% plotting dipoles
scatter3(x, y, z, 'filled');
axis('equal');
hold on;

end
