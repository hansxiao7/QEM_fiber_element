function f = constant_f(e, x,L,  N)

Sx = shape_fun(N, x, 1);
% [xn, wn] = lobatto_points(N);

B = Sx * e * 2/L;

f = (B(2)^2+1);