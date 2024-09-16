function max_eigenvectors = maxEigenvectors(KpcSym)
%MAXEIGENVECTORS

eigenvalues = eig(KpcSym);
max_eigenvalue = eigenvalues(end);

tol = 0.0001 * max_eigenvalue;
delta = max_eigenvalue - eigenvalues(end - 1);
i = 1;
while (delta < tol)
    i = i + 1;
    delta = max_eigenvalue - eigenvalues(end - i);
end

[max_eigenvectors, ~] = eigs(KpcSym, i, 'la');
end

