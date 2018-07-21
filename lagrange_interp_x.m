function Sx = lagrange_interp_x(N, z)
[x, ~] = lobatto_points(N);
for j = 1:N
    y = 0;
    n = length(x);
    for l=1:n
        if not(l==j)
            k = 1/(x(j)-x(l));
            for m=1:n
                if not(m==j) && not(m==l)
                    k = k*(z-x(m))/(x(j)-x(m));
                end
            end
            y = y + k;
        end
    end
    Sx(j) = y;
end

end