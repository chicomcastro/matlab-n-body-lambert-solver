%% script bla

clear;

% Do not change this (see showResults.m)
global shouldPlot;
shouldPlot = 1;
loadData;

%% Earth-Venus
% From 2B optim
earth_venus_transfer = fly(...
    108.9*24*60*60/ut, ... % Flight time Mars-Venus [ut]
    (R_earth_sun * [1 0 0] * Rz(0)), ...    % Departure from Earth [m]
    [-0.052512, 26.681, 0]... % Velocity to reach Mars orbit from Earth [km/s]
);
% TODO simular até cruzar a órbita de Vênus para refinar t_voo e x_t2

%% Venus swing by
% From 2B
% entro no swing_by com uma velocidade vinda da análise 2B
v_entra_swing_by = [-27.9215017466389,-24.7200847282275,0]; % [km/s]
% TODO pegar a velocidade do fim do estágio anterior

% To optimize
theta_soi_venus = 1.8102;   % [rad]
t_voo = 0.1;    % [ut]

% Simula o swing by
venus_pos = R_venus_sun * [1,0,0] * Rz(2.14795059794092);   % from 2B [m]
spaceship_pos = venus_pos + SOI_venus * [1,0,0] * Rz(theta_soi_venus); % [m]
venus_swing_by = swing_by(...
    t_voo, ...           % Simulation time (independant, it will stop when leaves SOI) [ut]
    spaceship_pos, ...   % Swing by start position [m]
    v_entra_swing_by,... % Velocity entering Venus SOI  [km/s]
    venus_pos...         % Venus position when spaceship entersits SOI [m]
);

error = (norm(venus_swing_by.error) - SOI_venus)/ud;
disp(error);

v_sai_swing_by = [-34.5525509338481,-21.7509589497666,0];
delta_v_swing_by = norm(v_sai_swing_by - venus_swing_by.v_t2/1000);

%% Venus-Mars
% From 2B optim
% tentando ver se saindo do swing by naturalmente alcança-se Marte
% utilizando as posições e velocidades calculadas no swing by
venus_mars_transfer = fly(...
    212.335011913553*24*60*60/ut, ... % Flight time Mars-Venus [ut]
    venus_swing_by.x_t2, ...    % Departure from Venus [m]
    venus_swing_by.v_t2/1000,... % Velocity to reach Mars orbit from Venus [km/s]
    venus_swing_by.target_x... % Venus phase after swing-by
);

%% Venus-Mars
% From 2B optim
% tenho que sair do swing-by com [-34.5525509338481,-21.7509589497666,0] de
% velocidade pra poder chegar em Marte
venus_mars_transfer = fly(...
    212.335011913553*24*60*60/ut, ... % Flight time Mars-Venus [ut]
    earth_venus_transfer.x_t2, ...    % Departure from Venus [m]
    [-34.5525509338481,-21.7509589497666,0],... % Velocity to reach Mars orbit from Venus [km/s]
    venus_swing_by.target_x... % Venus phase after swing-by
);

%%
shouldPlot = 0;
fun = @(x) perform_swing_by(x(1));
x0 = [pi/2];
options = optimset('MaxIter',1e4, 'MaxFunEvals', 1e4, 'TolFun', 1e-8, 'TolX', 1e-8);
x = fminsearch(fun,x0, options);
disp(x);
shouldPlot = 1;
%%
function result = fly(t_voo, initial_pos, initial_velocity_guess, venus_initial_pos)
    initial_state.t_voo = t_voo;                % [ut] (not [s], see normalization.m)
    initial_state.initial_pos = initial_pos;    % [m]
    initial_state.v0 = initial_velocity_guess;  % [km/s]

    x = initial_state;
    
    if nargin > 3
        x.venus_initial_pos = venus_initial_pos;
    end

    result = simulate(x);
end

function result = swing_by(t_voo, initial_pos, initial_velocity_guess, venus_initial_pos)    
    V_exit_earth = initial_velocity_guess;            % [km/s]
    
    if nargin <= 3
        loadData;
        theta_venus_f = 2.14795059794092;   % From 2B optim
        omega_venus_sun = norm([0,0,3.236706097374289e-07])*ut; % [rad/ut]
        theta_venus_i = theta_venus_f - omega_venus_sun * t_voo;
        venus_initial_pos = SOI_venus * [1,0,0] * Rz(theta_venus_i);
    end
    
    % Run all subroutines related to simulation
    main;
    
    global N;
    spaceship = y(:,3*(N-1)+1:3*(N-1)+3);
    planet = y(:,3*(N-2)+1:3*(N-2)+3);
    distance = spaceship - planet;
    distance_error = vecnorm(distance')' - SOI_venus/ud;
    for i=2:length(distance_error)
        if distance_error(i) > 0
            break;
        end
    end
    y = y(1:i-1,:);
    
    % The target will be my N-th body (state is based on N-1 position)
    r = y(end,3*(N-1)+1:3*(N-1)+3);     % [ud]
    v = y(end,3*(2*N-1)+1:3*(2*N-1)+3); % [ud/ut]
    
    % The target will be my (N-1)-th body (state is based on N-2 position)
    target = y(end,3*(N-2)+1:3*(N-2)+3);% [ud]
    error = (r(:)' - target(:)')*ud;    % [m] [row vector]
    
    if length(ie) > 0 && ie(end) == 1
        result.error = Inf;         % Do not count if ode has limited time
    else
        result.error = error;       % [m]
    end

    % Initial state
    result.x_t1 = y(1,3*(N-1)+1:3*(N-1)+3)*ud;      % [m]
    result.v_t1 = y(1,3*(2*N-1)+1:3*(2*N-1)+3)*ud/ut; % [m/s]
    
    % Final state
    result.x_t2 = r*ud;             % [m]
    result.v_t2 = v*ud/ut;          % [m/s]
    
    % Target state
    result.target_x = target*ud;      % [m]
    result.target_v = y(end,3*(2*N-2)+1:3*(2*N-2)+3)*ud/ut;      % [m/s]
end

function delta_v_swing_by = perform_swing_by(theta_entrada)
    loadData;
    t_voo = 0.1;
    v_entra_swing_by = [-27.9215017466389,-24.7200847282275,0]; % [km/s]

    % Simula o swing by
    initial_pos = R_venus_sun * [1,0,0] * Rz(2.14795059794092) + SOI_venus * [1,0,0] * Rz(theta_entrada);
    venus_swing_by = swing_by(...
        t_voo, ... % Flight time Mars-Venus [ut]
        initial_pos, ...    % Departure from Venus [m]
        v_entra_swing_by,... % Velocity to reach Mars orbit from Venus [km/s]
        2.14795059794092... % Venus phase before swing-by
    );

    v_sai_swing_by = [-34.5525509338481,-21.7509589497666,0];
    
    if t_voo <= 0
        delta_v_swing_by = Inf;
    else
        delta_v_swing_by = norm(v_sai_swing_by - venus_swing_by.v_t2/1000);
    end
end
