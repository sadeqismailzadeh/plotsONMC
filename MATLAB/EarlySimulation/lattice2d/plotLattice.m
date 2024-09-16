function plotLattice(R, r, l)
% plot a visualization of lattice including cell boundaries
% 
% R           3 * 2 matrix of lattice vectors
% r           3 * n  array of relettive particle positions inside a cell
% l           lattice size

R1 = R(:,1);
R2 = R(:,2);

N = l*l;
n = size(r,2);  %number of dipoles inside a cell

x = zeros(1, n*N);
y = zeros(1, n*N);
z = zeros(1, n*N);

% storing location of dipoles in x and y arrays
a = 1;
for i = 0:l-1
    for j = 0:l-1
        for k = 1: n
            x(a) = i*R1(1)+j*R2(1)+r(1,k);
            y(a) = i*R1(2)+j*R2(2)+r(2,k);
            z(a) = r(3, k);
            
            a = a + 1;
        end
    end
end

% plotting dipoles
scatter3(x, y, z,'filled');
axis('equal');
hold on;

% plotting cell boundaries
for i = 0:l
    line([i*R2(1)  i*R2(1)+l*R1(1)],[i*R2(2)  i*R2(2)+l*R1(2)],'Color','r')
end

for i = 0:l
    line([i*R1(1)  i*R1(1)+l*R2(1)],[i*R1(2)  i*R1(2)+l*R2(2)],'Color','r')
end

end
