function plotCell(R ,r)
% plot a visualization of a unit cell
% r      3 * n  matrix of relettive particle positions inside a cell

n = size(r,2);  %number of dipoles inside a cell

R1 = R(:,1);
R2 = R(:,2);
R3 = R(:,3);

x = zeros(1, n);
y = zeros(1, n);
z = zeros(1, n);

a = 1;

for t = 1: n
    x(a) = r(1,t);
    y(a) = r(2,t);
    z(a) = r(3,t);
    
    a = a + 1;
end

scatter3(x, y, z, 'filled');
% xlim([0, R1(1)+R2(1)+R3(1)])
% ylim([0, R1(2)+R2(2)+R3(2)])
% zlim([0, R1(3)+R2(3)+R3(3)])
axis('equal');

end