N = 5;

xn = sym('x',N);
wn = sym('w',N);

syms x;
syms L;
syms rho;
syms A;

for i = 1:N
    C(i) = sym(1);
    for j = 1:N
        if i~=j
            C(i) = C(i) * (x - xn(j))/(xn(i) - xn(j));
        end
    end
end

Cx = diff(C, x);

hi0 = sym(zeros(1, N));

h11 = (x^2 - (xn(1) + xn(N))*x + xn(1)*xn(N))/(xn(1) - xn(N)) * C(1);
hn1 = (x^2 - (xn(1) + xn(N))*x + xn(1)*xn(N))/(xn(N) - xn(1)) * C(N);
h10 = (x - xn(N)) / (xn(1) - xn(N)) * (1- (x-xn(1))*(subs(Cx(1), x, xn(1))+ 1/ (xn(1) - xn(N)))) * C(1);
hn0 = (x - xn(1)) / (xn(N) - xn(1)) * (1- (x-xn(N))*(subs(Cx(N), x, xn(N))+ 1/ (xn(N) - xn(1)))) * C(N);

for i = 2:N-1
    hi0(i) = (x - xn(1)) * (x - xn(N))/ ((xn(i) - xn(1)) * (xn(i) - xn(N))) * C(i);
end
hi0(1) = h10;
hi0(N) = hn0;

S = sym(zeros(2, 2*N + 2));
for i = 1:N
    S(1, 2*i) = C(i);
end

S(2,1) = h11;
S(2, 2*N +2) = hn1;

for j = 1:N
    S(2, 2*j+1) = hi0(j);
end
S(1,1) = S(1,2);
S(1,2) = 0;
S(2,2) = S(2,3);
S(2,3) = S(2,1);
S(2,1) = 0;
M = simplify(int(rho* A * L / 2 * (S') * S, -1, 1));

matlabFunction(M, 'file', 'consistent_mass_5.m');
matlabFunction(S, 'file', 'symbol_S.m');