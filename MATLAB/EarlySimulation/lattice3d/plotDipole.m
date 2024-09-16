function plotDipole(R, A, L)
% Plots a dipole in position R i.e. an arrow whose center is at R
% tip of arrow is red
% 
% input:
% R      position of dipole   
% A      spherical Angles (theta = A(1), phi = A(2))
% L      length (indicates strength of dipole)

x = zeros(3,1);
y = zeros(3,1);

th = A(1);
phi = A(2);

x0 = R(1);
y0 = R(2);
z0 = R(3);

xi = x0 - L*sin(th)*cos(phi);
yi = y0 - L*sin(th)*sin(phi);
zi = z0 - L*cos(th);

xm = x0 + (L-L*0.3)*sin(th)*cos(phi);
ym = y0 + (L-L*0.3)*sin(th)*sin(phi);
zm = z0 + (L-L*0.3)*cos(th);

xf = x0 + L*sin(th)*cos(phi);
yf = y0 + L*sin(th)*sin(phi);
zf = z0 + L*cos(th);


line([xi  xm],[yi  ym],[zi  zm], 'linewidth', 1, 'Color', 'b')
line([xm  xf],[ym  yf],[zm  zf], 'linewidth', 4, 'Color', 'r')

end

