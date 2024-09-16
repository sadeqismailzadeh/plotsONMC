ls = [1 2 3 4 5 6 7 8];
invls = 1 ./ ls;
vals = [4.516 3.928 3.606 3.399 3.254 3.1455 3.0607 2.9924];
% plot(invls, vals)
 interp1(invls, vals,  0, 'spline',  'extrap')