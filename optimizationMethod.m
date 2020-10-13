%% Inputs
% Define optimization parameters
max_iteration = 1000;                                  % Iteractions to run
num_particles = 1000;                                  % Particles in swarm

clear lower_boundary upper_boundary

% Do not change this (see showResults.m)
global shouldPlot;
shouldPlot = 0;

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
delta_v = [ 0 0.5 0];
lower_boundary.v0 = [36.5891 55.5695 0]*0;
upper_boundary.v0 = [36.5891 55.5695 0]*2;

%%
lower_boundary = struct_2_state(lower_boundary); % Lower boundary
upper_boundary = struct_2_state(upper_boundary); % Upper boundary

%% Optimization
% Run optimization
tentativa = 0;
while tentativa < 1
    disp("Tentativa: " + tentativa);
    disp("Rodando para num_particles = " + num_particles);
    tic
    pso;
    execution_time = toc;
    disp("Elapsed time was " + execution_time + " seconds.")

    x = state_2_struct(best_global);
    if x.t_voo ~= 0
        result = simulate(x);
        disp("Minimum error was " + norm(result.error)/ud + " [ud]");
    else
        disp("No non-zero best_global evaluated");
    end
    
    tentativa = tentativa + 1;
end
    
function x = struct_2_state(struct)
    x = [
        struct.t_voo, ...
        struct.initial_pos, ...
        struct.target_pos, ...
        struct.v0 ...
    ];
end

function x = state_2_struct(state)
    x.t_voo = state(1);
    x.initial_pos = state(2:4);
    x.target_pos = state(5:7);
    x.v0 = state(8:10);
end
