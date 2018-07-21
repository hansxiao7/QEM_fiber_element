function [MM, Qint,new_sections] = internal_force_inelastic(e, L, N, sections)

[xn, wn] = lobatto_points(N);

Qint = zeros(2*N+2, 1);


for i = 1:N
    x = xn(i);
    w = wn(i);
    Sx = shape_fun(N, x, 1)*2/L;
    Sxx = shape_fun(N,x, 2) * 4/L^2;
    strain_x = axial_strain(e, x, L, N);
    if isreal(strain_x) ~= 1
        error('strain_x negative');
    end
    phi = curvature(e, x, L, N);
    [moment, axial_force, new_sections(i)] = section_analysis(sections(i), phi, strain_x);
    if i == 1
        MM = moment;
    end
    Qint = Qint + w * L / 2 * (moment * Sxx(2,:)' + axial_force * (Sx(1,:)' + Sx(1,:)'*Sx(1,:)*e + Sx(2,:)'*Sx(2,:)*e));
end