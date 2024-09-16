function plotLatticeDipole(R, r, L, A, l)
% PLOTLATTICEDIPOLE plots a visualization of lattice with dipoles 
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

% plot dipoles
for i = 0:L-1
    for j = 0:L-1
        for k = 0:L-1
            
            for t = 1: n
                x = i*R1(1) + j*R2(1) + k*R3(1)+ r(1,t);
                y = i*R1(2) + j*R2(2) + k*R3(2)+ r(2,t);
                z = i*R1(3) + j*R2(3) + k*R3(3)+ r(3,t);
                
                plotDipole([x, y, z], A(t, :), l(t));
                
            end
            
        end
    end
end

axis equal
% hold on;

end