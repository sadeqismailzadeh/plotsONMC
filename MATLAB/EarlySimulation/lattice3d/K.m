function K_out = K(R_cells, r ,n)
% interaction of two cells with the relative distance of R_cells
% 
% R_cells     distance between two cells  
% n           number of dipoles inside a cell
% r           3 * n  array of relettive particle positions inside a cell

K_out = zeros(n, n);  
    for i = 1:n
        for j = 1:n
            K_out(3*i-2:3*i, 3*j-2:3*j) = J(-r(:,i) + R_cells + r(:,j));
        end
    end
end
