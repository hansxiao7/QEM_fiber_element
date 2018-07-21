function section = sections_init(N)
%this function is to define the element sections. Each element has same
%section properties, and have N sections (Lobatto points) for analysis

%n_sec = N;
for m = 1:N
        section(m).n_conc_core = 10;
        section(m).n_conc_cover = 12;
        section(m).n_conc = section(m).n_conc_core + section(m).n_conc_cover; %number of concrete fibers
        section(m).n_steel = 2;
        section(m).y_conc_core = [9; 7; 5; 3; 1; -1; -3; -5; -7; -9];
        section(m).y_conc_cover = [11.5; 9; 7; 5; 3; 1; -1; -3; -5; -7; -9; -11.5];
        section(m).y_conc = [section(m).y_conc_core; section(m).y_conc_cover];
        section(m).y_steel = [10; -10];

        section(m).A_conc_core = [28*2; 28*2; 28*2; 28*2; 28*2; 28*2; 28*2; 28*2; 28*2; 28*2];
        section(m).A_conc_cover = [34*3; 6*2; 6*2; 6*2;6*2;6*2;6*2;6*2;6*2;6*2;6*2;34*3];
        section(m).A_conc = [section(m).A_conc_core; section(m).A_conc_cover];
        section(m).A_steel = [3.14; 3.14];

        %fiber concrete properties
        fc_cover = -5.2;
        strain_0_cover = -0.002;
        strain_u_cover = -0.018;
        Z_cover = 420;
        
        fc_core = -5.6;
        strain_0_core = -0.0022;
        strain_u_core = -0.018;
        Z_core = 53.485;

        %fiber reinf properties-- built in steel01 model


        %define the fibers
        %concrete
        for j = 1:section(m).n_conc
            if j <= section(m).n_conc_core
                section(m).conc_fiber(j).fc = fc_core;
                section(m).conc_fiber(j).strain_0 = strain_0_core;
                section(m).conc_fiber(j).strain_u = strain_u_core;
                section(m).conc_fiber(j).Z = Z_core;
            else
                section(m).conc_fiber(j).fc = fc_cover;
                section(m).conc_fiber(j).strain_0 = strain_0_cover;
                section(m).conc_fiber(j).strain_u = strain_u_cover;
                section(m).conc_fiber(j).Z = Z_cover;
            end
            section(m).conc_fiber(j).strain_history = [0,0];
            section(m).conc_fiber(j).stress_history = [0,0];
            section(m).conc_fiber(j).strain_r = 0;
            section(m).conc_fiber(j).strain_p = 0;
            section(m).conc_fiber(j).E = 0;
        end
        
        %reinf
        for k = 1:section(m).n_steel
            section(m).steel_fiber(k).strain_history = [0,0];
            section(m).steel_fiber(k).stress_history = [0,0];
            section(m).steel_fiber(k).ita = 0;
            section(m).steel_fiber(k).E = 0;
        end

end

%initialize the tangent modulus
for i = 1:N
    [~, ~, section(i)] = section_analysis(section(i), 0, 0);
end

end