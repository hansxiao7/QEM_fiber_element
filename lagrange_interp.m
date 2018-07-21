function S = lagrange_interp(N, x)
[xn, ~] = lobatto_points(N);
for i = 1:N
    S(i) = 1;
    for j = 1:N
        if i~=j
            S(i) = S(i) * (x - xn(j))/(xn(i) - xn(j));
        end
    end
end


end