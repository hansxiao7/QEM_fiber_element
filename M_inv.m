function MM = M_inv(rho, A, L, N)

M = mass(rho, A, L, N);

MM = zeros(2*N+2, 2*N+2);

for i = 1:2*N+2
    if M(i, i) ~= 0
        MM(i, i) = 1 / M(i,i);
    end
end

end