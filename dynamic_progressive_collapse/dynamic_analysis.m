clc
clear

load('gravity_result_7.mat');
Qg = gravity_force(rho_g, N, L, A);

%progressive collapse beam condition
Qg(1:3,1) = zeros(3,1);
Qg(2*N, 1) = 0;
Qg(2*N+2, 1) = 0;

%external force from column
Q_column = zeros(2*N+2, 1);
Q_column(2*N + 1, 1) = -4.95;

%total external force = beam gravity + column gravity
Qext = Qg + Q_column;

% Mass matrix (constant)
rho = 2.17e-7;
M = mass(rho, A, L, N)+rotation_mass(rho, I, L, N);
%additional mass from column
M(2*N+1, 2*N+1) = M(2*N+1, 2*N+1) + 0.0128;
M(2*N, 2*N) = M(2*N, 2*N) + 0.0128;

Minv = inv(M);
%left end fiexed
Minv(1:3, :) = zeros(3, 2*N+2);
Minv(:, 1:3) = zeros(2*N+2, 3);
%right end rotation and horizontal fixed
Minv(2*N, :) = zeros(1, 2*N+2);
Minv(2*N+2, :) = zeros(1, 2*N+2);
Minv(:, 2*N) = zeros(2*N + 2, 1);
Minv(:, 2*N+2) = zeros(2*N + 2, 1);

%deicde the delta_t
delta_t = 8.65e-06;
tend = 0.5;

%damping
alpha = 0.01;

%set the saving space for analysis
nout = floor(tend/delta_t) + 1;
data.t = zeros(1, nout);
data.e = zeros(2*N + 2, nout);
data.Q = zeros(2*N + 2, nout);
data.M = zeros(1, nout);
data.strain_ends = zeros(2, nout);
data.curvature_ends = zeros(2, nout);
data.reinf_stress = zeros(4, nout);
data.reinf_strain = zeros(4, nout);
data.conc_stress = zeros(4, nout);

data.t(1) = 0;
data.e(:,1) = e;

old_e = e;

% Loop over all output times.
hw = waitbar(0,'Initializing waitbar...');

for j = 1:nout-1
    [MM, Qint,sections] = internal_force_inelastic(e, L, N, sections);
    data.Q(:,j) = Qint; 
    data.M(j) = MM;
    data.strain_ends(1,j) = axial_strain(e, -1, L, N);
    data.strain_ends(2,j) = axial_strain(e, 1, L, N);
    data.curvature_ends(1,j) = curvature(e, -1, L, N);
    data.curvature_ends(2,j) = curvature(e, 1, L, N);
    data.reinf_stress(1,j) = sections(1).steel_fiber(1).stress_history(2);
    data.reinf_stress(2,j) = sections(1).steel_fiber(2).stress_history(2);
    data.reinf_stress(3,j) = sections(N).steel_fiber(1).stress_history(2);
    data.reinf_stress(4,j) = sections(N).steel_fiber(2).stress_history(2);
    
    data.reinf_strain(1,j) = sections(1).steel_fiber(1).strain_history(2);
    data.reinf_strain(2,j) = sections(1).steel_fiber(2).strain_history(2);
    data.reinf_strain(3,j) = sections(N).steel_fiber(1).strain_history(2);
    data.reinf_strain(4,j) = sections(N).steel_fiber(2).strain_history(2);
    
    data.conc_stress(1,j) = sections(1).conc_fiber(11).stress_history(2);
    data.conc_stress(2,j) = sections(1).conc_fiber(22).stress_history(2);
    data.conc_stress(3,j) = sections(N).conc_fiber(11).stress_history(2);
    data.conc_stress(4,j) = sections(N).conc_fiber(22).stress_history(2);
    
    %boundary condition
    Qint(1:3, 1) = zeros(3,1);
    Qint(2*N, 1) = 0;
    Qint(2*N + 2, 1) = 0;
    
    temp = e;
    e = delta_t^2 * Minv * (Qext - Qint) +2*e - old_e - delta_t* alpha * (e - old_e);
    e(1:3,1) = zeros(3,1);
    e(2*N, 1) = 0;
    e(2*N + 2, 1) = 0;
    old_e = temp;
    data.e(:,j+1) = e;
    data.t(j+1) = data.t(j) + delta_t;
    t_str = sprintf('t = %.4f', j*delta_t);
    waitbar(j/(nout-1),hw,t_str);
end
close(hw);

