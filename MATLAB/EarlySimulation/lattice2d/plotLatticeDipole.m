function plotLatticeDipole(R, r, L, A, l)
% PLOTLATTICEDIPOLE  plots a visualization of lattice with dipoles including cell boundaries
% 
% R          3 * 2 matrix of lattice vectors
% r          3 * n  array of relettive particle positions inside a cell
% L          lattice size
% A          array of spherical Angles of dipoles,
%  			 theta of ith  dipole = A(1,i), phi of ith dipole = A(2,i)
% l          lengths of dipoles

R1 = R(:,1);
R2 = R(:,2);

n = size(r,2);  %number of dipoles inside a cell

% plot dipoles
a = 1;
for i = 0:L-1
    for j = 0:L-1
        for k = 1: n
            x = i*R1(1)+j*R2(1)+r(1,k);
            y = i*R1(2)+j*R2(2)+r(2,k);
            z = r(3, k);
            
            plotDipole([x, y, z], A(k,:), l(k));
        end
    end
end


axis('equal');
hold on;

% plot cell boundaries
for i = 0:L
    line([i*R2(1)  i*R2(1)+L*R1(1)],[i*R2(2)  i*R2(2)+L*R1(2)],'Color','r')
end

for i = 0:L
    line([i*R1(1)  i*R1(1)+L*R2(1)],[i*R1(2)  i*R1(2)+L*R2(2)],'Color','r')
end
end