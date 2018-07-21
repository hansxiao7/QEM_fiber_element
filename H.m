function [h11, hn1, hi0] = H(N, x)

S = lagrange_interp(N, x);

[xn, ~] = lobatto_points(N);

Sx1 = lagrange_interp_x(N, xn(1));
Sx1 = Sx1(1);
SxN = lagrange_interp_x(N, xn(N));
SxN = SxN(N);

h11 = (x^2 - (xn(1) + xn(N))*x + xn(1) * xn(N)) / (xn(1) - xn(N)) * S(1);
hn1 = (x^2 - (xn(1) + xn(N))*x + xn(1) * xn(N)) / (xn(N) - xn(1)) * S(N);

h10 = (x - xn(N)) / (xn(1) - xn(N)) * (1 - (x - xn(1)) * (Sx1 + 1 / (xn(1) - xn(N)))) * S(1);
hN0 = (x - xn(1)) / (xn(N) - xn(1)) * (1 - (x - xn(N)) * (SxN + 1 / (xn(N) - xn(1)))) * S(N);

hi0 = zeros(1,N);

for i = 2:N-1
    hi0(i) = (x - xn(1)) * (x - xn(N)) / ((xn(i) - xn(1)) * (xn(i) - xn(N))) * S(i);
end

hi0(1) = h10;
hi0(N) = hN0;

end
