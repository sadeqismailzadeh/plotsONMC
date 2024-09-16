function [J_out] = J(r)
% interaction of two dipoles with the relative distance of r

x = r(1);
y = r(2);
z = r(3);
J_out = [2*x^2-y^2-z^2      3*x*y              3*x*z;
         3*x*y              2*y^2-x^2-z^2      3*y*z;
         3*x*z              3*y*z              2*z^2-x^2-y^2] / ((x^2 + y^2 + z^2)^2.5);
     
end