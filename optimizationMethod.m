%% Inputs
% Define optimization parameters
max_iteration = 100;                                  % Iteractions to run
num_particles = 100;                                  % Particles in swarm

clear lower_boundary upper_boundary

%% Pré declarations and calculations
% This is the section that can be security changed to fit your problem
% Others section can cause problems to change, then it's your risk
loadData;

t_voo = 4.8200;   % Simulation time [rad]
initial_pos = (R_earth_sun * [1 0 0] + altitude_from_earth * [1 0 0]);	% [m]
target_pos = R_mars_sun * [-1 0 0];     % [m]

%% Boundaries struct

% Fixed values
lower_boundary.t_voo = t_voo;
upper_boundary.t_voo = t_voo;
lower_boundary.initial_pos = initial_pos;
upper_boundary.initial_pos = initial_pos;
lower_boundary.target_pos = target_pos;
upper_boundary.target_pos = target_pos;

% Variable values
lower_boundary.v0 = [0, V_oe_earth + V_earth_sun / 2, 0];
upper_boundary.v0 = ones(1,3)*V_oe_inertial;

%%
lower_boundary = struct_2_boundary(lower_boundary); % Lower boundary
upper_boundary = struct_2_boundary(upper_boundary); % Upper boundary

%% Optimization
% Run optimization
tentativa = 0;
while tentativa < 1
    disp("Tentativa: " + tentativa);
    disp("Rodando para num_particles = " + num_particles);
    tic
    pso;
    execution_time = toc;
    disp("Elapsed time was " + execution_time + "seconds.")
    tentativa = tentativa + 1;
end
    
function x = struct_2_boundary(struct)
    x = [
        struct.t_voo, ...
        struct.initial_pos, ...
        struct.target_pos, ...
        struct.v0 ...
    ];
end
