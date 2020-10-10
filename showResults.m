% function showResults()
% Plot trajectory and observables quantities after running any of the
% repository N-body simulation / optimization

global shouldPlot;
shouldPlot = 1;

% Define custom initial conditions
t_voo = x.t_voo;                % [ut]
initial_pos = x.initial_pos;    % [m]
V_exit_earth = x.v0;            % [km/s]

% Call main routine to call subroutines
main;
