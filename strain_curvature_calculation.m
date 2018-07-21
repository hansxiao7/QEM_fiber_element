strain_x = [];
xx= [];
phi = [];
for i = 1:2001
    x = -1 + (i-1)*0.001;
    xx = [xx;x];
    strain_x = [strain_x; axial_strain(e, x, L, N)];
    phi = [phi; curvature(e,x,L,N)];
end
strain_x = [];
xx= [];
phi = [];
for j = 1:13
    x = xn(j);
    xx =[xx;x];
    strain_x = [strain_x; axial_strain(e, x, L, N)];
    phi = [phi; curvature(e,x,L,N)];
end