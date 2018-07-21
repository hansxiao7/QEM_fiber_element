xx = [];
vv = [];
for i = 1:201
    x = -1 + (i-1)*0.01;
    xx = [xx;x];
    S = shape_fun(N, x, 0);
    v = S * e;
    vv= [vv; v(1)];
end
xx = xx * L/2 + L/2;