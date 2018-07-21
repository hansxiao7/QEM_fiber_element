clc
clear
x= [];
y=[];
E = [];
result = [];
ita_result = [];
max = [];
min = [];
strain_0 = [];
for i = 1:674
    if i == 1
        steel.strain_history = [0,0];
        steel.stress_history = [0,0];
        steel.strain_r = -69/29000;
        steel.max_strain_r = 69/29000;
        steel.min_strain_r = -69/29000;
        steel.stress_r = -69;
        steel.strain_0 = 69/29000;
        steel.stress_0 = 69;
        steel.ita = 0;
    end
    steel = steel01(strain(i), steel.strain_history, steel.stress_history, steel.ita...
    , steel.strain_r, steel.stress_r, steel.strain_0, steel.stress_0, steel.max_strain_r, steel.min_strain_r);
%     x = [x, steel.strain_history(2)];
    strain_0 = [strain_0, steel.strain_0];
    E = [E, steel.E];
    result = [result, steel.stress];
    ita_result = [ita_result, steel.ita];
    max = [max, steel.max_strain_r];
    min = [min, steel.min_strain_r];
end

for j = 1:10
    steel = steel01(0.0005-0.0001 * j, steel.strain_history, steel.stress_history, steel.ita...
    , steel.strain_r, steel.stress_r, steel.strain_0, steel.stress_0, steel.max_strain_r, steel.min_strain_r);
    x = [x, steel.strain_history(2)];
    strain_0 = [strain_0, steel.strain_0];
    result = [result, steel.stress];
    E = [E, steel.E];
    ita_result = [ita_result, steel.ita];
     max = [max, steel.max_strain_r];
    min = [min, steel.min_strain_r];
end
steel.stress
for k = 1:10
        steel = steel01(-0.0005+0.0001 * k, steel.strain_history, steel.stress_history, steel.ita...
    , steel.strain_r, steel.stress_r, steel.strain_0, steel.stress_0, steel.max_strain_r, steel.min_strain_r);
    x = [x, steel.strain_history(2)];
    result = [result, steel.stress];
    E = [E, steel.E];
    ita_result = [ita_result, steel.ita];
     max = [max, steel.max_strain_r];
    min = [min, steel.min_strain_r];
end
steel.stress
y = result;