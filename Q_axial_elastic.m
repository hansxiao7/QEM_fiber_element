function Fint = Q_axial_elastic(e, x, L, N, E, A)
% -1 <= x <= 1

Sx = shape_fun(N, x, 1);

Sx1 = Sx(1,:);
Sx2 = Sx(2,:);

Fint = E*A*L/2 * (Sx1'*Sx1 * e * 4/L^2);

end