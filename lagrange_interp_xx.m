function Sxx = lagrange_interp_xx(N, x)
[xn, ~] = lobatto_points(N);
for i = 1:N
    Sxx(i) = 0;
    for l = 1:N
        y = 0;
        if l~=i
            for m = 1:N
                if m~=i && m~=l
                    temp = 1/(xn(i)-xn(m));
                    for k = 1:N
                        if k~=i && k~=l && k~=m
                            temp = temp * (x-xn(k))/(xn(i)-xn(k));
                        end
                    end
                    y = y + temp;
                end
            end
            Sxx(i) = Sxx(i) + y / (xn(i)-xn(l));
        end
    end
end

end

                    