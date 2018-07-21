clc
clear

load('dynamic_results.mat');
load('gravity_result.mat','sections');
stress = [];
strain = [];
times = j;
for i = 1:times
    i
    e = data.e(:,i);
    stress = [stress;sections(1).steel_fiber(1).stress_history(2)];
    strain = [strain;sections(1).steel_fiber(1).strain_history(2)];
    [Qint,sections] = internal_force_inelastic(e, L, N, sections);
end
    