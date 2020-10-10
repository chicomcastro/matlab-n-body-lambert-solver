clear;
global shouldPlot;
shouldPlot = 0;

%%
loadData;
V_oe_earth = sqrt(G_metric*M_earth/altitude_from_earth);     % [m/s]
V_oe_inertial = V_oe_earth + V_earth_sun;
V_oe_inertial = V_oe_inertial/1e3;  % to [km/s]

t_voo = 4.8200;   % Tempo [rad]
semieixo = (R_mars_sun + R_earth_sun) / 2;
V_exit_earth = sqrt(G_metric*M_sun*(1/semieixo+2/R_earth_sun))/1000*[0 1 0]; %km/s
initial_velocity_guess = [36.5891 55.5695 0];  % km/s

initial_pos = (R_earth_sun * [1 0 0] + altitude_from_earth * [1 0 0]);  % m
target_pos = R_mars_sun * [-1 0 0]; % m

%% Initial state struct
initial_state.t_voo = t_voo;
initial_state.initial_pos = initial_pos;
initial_state.target_pos = target_pos;
initial_state.v0 = initial_velocity_guess;

x = initial_state;
%%
itr = 1;
delta = 1e-2;
error = [];
erro_dist = ones(1,3)*ud;
disp("Iniciando busca de uma solução...");
tic
while norm(erro_dist)/ud > 5e-2
    if mod(itr, 10) == 0
        disp("---");
        disp("Iteração: " + itr);
        disp("Error: " + norm(erro_dist)/ud);
    end
    
    itr = itr + 1;
    x1 = x;
    result = coust(x1);
    
    % vx variation
    xx = x;
    xx.v0 = xx.v0 + delta*[1,0,0];
    result_dx = coust(xx);
    dxdvx = (result_dx.x_t2 - result.x_t2)/delta;
    
    % vy variation
    xy = x;
    xy.v0 = xy.v0 + delta*[0,1,0];
    result_dy = coust(xy);
    dxdvy = (result_dy.x_t2 - result.x_t2)/delta;
    
    % vz variation
    xz = x;
    xz.v0 = xz.v0 + delta*[0,0,1];
    result_dz = coust(xz);
    dxdvz = (result_dz.x_t2 - result.x_t2)/delta;
    
    % variation matrix
    dXdV = [dxdvx; dxdvy; dxdvz]';                  % [s]
    erro_dist = result.error;                       % [m]
    delta_velocidade = -1*dXdV\erro_dist(:);        % [m/s]
    x.v0 = x.v0 + delta_velocidade(:)'/1000;        % [km/s]
end
toc
disp("Solução encontrada!");
disp(x);
disp("Custo de saída: " + abs(norm(x.v0) - V_oe_inertial));
