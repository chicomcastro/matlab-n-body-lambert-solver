%% script bla
% inputs:
% V_sai_terra = [-5.1129, 27.293, 0] [km/s]
% theta_venus = 3.1405 [rad]
% theta_soi_venus [rad]
% fracao_impulso
% magnitude_impulso [km/s]

% Do not change this (see showResults.m)
global shouldPlot;
shouldPlot = 0;
loadData;

%% Earth-Venus
% From 2B optim
V_sai_terra = [-5.1129       27.293            0];    % [km/s]
earth_venus_transfer = fly(...
    121.12*24*60*60/ut, ... % Flight time Mars-Venus [ut]
    (R_earth_sun * [1 0 0] * Rz(0)), ...    % Departure from Earth [m]
    V_sai_terra... % Velocity to reach Mars orbit from Earth [km/s]
);
V_inicial = V_earth_sun * [0,1,0]/1000;    % [km/s]
delta_v_saida_terra = norm(V_sai_terra - V_inicial);    % [km/s]

v_inf = norm(V_sai_terra*1000 - V_earth_sun*[0,1,0]);   % V_saida em relação à Terra[m/s]
v_saida = sqrt(norm(v_inf)^2 + 2*G_metric*M_earth/(altitude_from_earth + R_earth)); % V_p que vai dar V_inf na saída hiperbólica
delta_v_saida_terra_2B = norm(v_saida - V_oe_earth)/1000;   % delta V de saída por 2B [km/s]

% TODO simular até cruzar a órbita de Vênus para refinar t_voo e x_t2

%% Venus swing by
% From 2B
% entro no swing_by com uma velocidade vinda da análise 2B
%v_entra_swing_by = [-5.1471      -37.727            0]; % [km/s]
v_entra_swing_by = earth_venus_transfer.v_t2/1000; % [km/s]
% TODO pegar a velocidade do fim do estágio anterior

% Can be anything sufficiently big to finish 1 non-propulsed swing by
t_voo = 0.1;    % [ut]

% To optimize
theta_soi_venus = 0.52135;   % [rad]
if ~exist('fracao_impulso', 'var')
    fracao_impulso = 0.512;
end
if ~exist('magnitude_impulso', 'var')
    magnitude_impulso = 0.85;
end

% Simula o swing by
% Avaliação do swing by
venus_pos = R_venus_sun * [1,0,0] * Rz(3.1405);   % from 2B [m]
spaceship_pos = venus_pos + SOI_venus * [1,0,0] * Rz(theta_soi_venus); % [m]
venus_swing_by = swing_by(...
    t_voo, ...           % Simulation time (independant, it will stop when leaves SOI) [ut]
    spaceship_pos, ...   % Swing by start position [m]
    v_entra_swing_by,... % Velocity entering Venus SOI  [km/s]
    venus_pos...         % Venus position when spaceship entersits SOI [m]
);
t_swing_by = venus_swing_by.t(end); % [ut]

% Pré-impulso
venus_swing_by_pre = swing_by(...
    t_swing_by*fracao_impulso, ...           % Simulation time (independant, it will stop when leaves SOI) [ut]
    spaceship_pos, ...   % Swing by start position [m]
    v_entra_swing_by,... % Velocity entering Venus SOI  [km/s]
    venus_pos...         % Venus position when spaceship entersits SOI [m]
);

% Pós-impulso
v_pre_impulso = venus_swing_by_pre.v_t2/1000;
v_pos_impulso = v_pre_impulso + magnitude_impulso * v_pre_impulso/norm(v_pre_impulso);
venus_swing_by_pos = swing_by(...
    t_swing_by*(1-fracao_impulso), ...           % Simulation time (independant, it will stop when leaves SOI) [ut]
    venus_swing_by_pre.x_t2, ...   % Swing by start position [m]
    v_pos_impulso,...    % Velocity entering Venus SOI  [km/s]
    venus_swing_by_pre.target_x...         % Venus position when spaceship entersits SOI [m]
);

%v_sai_swing_by = [0.51576      -40.787            0];    % [km/s]
%v_sai_swing_by = venus_swing_by.v_t2/1000;    % [km/s]
%delta_v_swing_by = norm(v_sai_swing_by - venus_swing_by_pos.v_t2/1000);    % [km/s]
delta_v_swing_by = abs(magnitude_impulso);    % [km/s]

%% Venus-Mars
% From 2B optim
% tentando ver se saindo do swing by naturalmente alcança-se Marte
% utilizando as posições e velocidades calculadas no swing by
venus_mars_transfer = fly(...
    212.335011913553*24*60*60/ut, ... % Flight time Mars-Venus [ut]
    venus_swing_by_pos.x_t2, ...    % Departure from Venus (SC position after swing by) [m]
    venus_swing_by_pos.v_t2/1000,... % Velocity to reach Mars orbit from Venus [km/s]
    venus_swing_by_pos.target_x... % Venus position after swing-by
);

%% Venus-Mars
% From 2B optim
% tenho que sair do swing-by com [0.51576      -40.787            0] de
% velocidade pra poder chegar em Marte
%venus_mars_transfer = fly(...
%    212.335011913553*24*60*60/ut, ... % Flight time Mars-Venus [ut]
%    venus_swing_by.x_t2, ...    % Departure from Venus [m]
%    [0.51576      -40.787            0],... % Velocity to reach Mars orbit from Venus [km/s]
%    venus_swing_by.target_x... % Venus phase after swing-by
%);

%%
tangent_direction = cross([0,0,1], venus_mars_transfer.x_t2(:)')/norm(cross([0,0,1], venus_mars_transfer.x_t2(:)'));
V_final = V_mars_sun * tangent_direction;    % [m/s]
delta_v_chegada_marte = norm(venus_mars_transfer.v_t2 - V_final)/1000;    % [km/s]

v_inf = norm(venus_mars_transfer.v_t2) - V_mars_sun;
v_chegada = sqrt(norm(v_inf)^2 + 2*G_metric*M_mars/(altitude_from_mars + R_mars));
delta_v_chegada_marte_2B = norm(V_oe_mars - v_chegada)/1000;   % delta V de saída por 2B [km/s]

if venus_mars_transfer.error > 0.01
    delta_v_chegada_marte = NaN;
    delta_v_chegada_marte_2B = NaN; 
end

%%
delta_v = [...
    delta_v_saida_terra;...
    delta_v_swing_by;...
    delta_v_chegada_marte...
];
delta_v_2B = [...
    delta_v_saida_terra_2B;...
    delta_v_swing_by;...
    delta_v_chegada_marte_2B...
];
delta_v_total = sum(delta_v);
delta_v_total_2B = sum(delta_v_2B);

%%
function result = fly(t_voo, initial_pos, initial_velocity_guess, venus_initial_pos)    
    V_exit_earth = initial_velocity_guess;            % [km/s]
    
    if nargin <= 3
        loadData;
        theta_venus_f = 3.1405;   % From 2B optim
        omega_venus_sun = norm([0,0,3.236706097374289e-07])*ut; % [rad/ut]
        theta_venus_i = theta_venus_f - omega_venus_sun * t_voo;
        venus_initial_pos = R_venus_sun * [1,0,0] * Rz(theta_venus_i);
    end
    
    % Run all subroutines related to simulation
    main;
    
    global N;
    % Verify it have crossed Mars' orbit
    spaceship_pos = y(:,3*(N-1)+1:3*(N-1)+3);
    distance_error = vecnorm(spaceship_pos')' - R_mars_sun/ud;
    for i=1:length(distance_error)
        if distance_error(i) > 0
            break;
        end
    end
    if i < length(distance_error)
        y = y(1:i-1,:);
        t = t(1:i-1);
    end
    
    % The spaceship will be my N-th body (state is based on N-1 position)
    r = y(end,3*(N-1)+1:3*(N-1)+3);     % [ud]
    v = y(end,3*(2*N-1)+1:3*(2*N-1)+3); % [ud/ut]
    
    % The target will be my (N-1)-th body (state is based on N-2 position)
    target = y(end,3*(N-2)+1:3*(N-2)+3);% [ud]
    
    if (length(ie) > 0 && ie(end) == 1)
        result.error = Inf;         % Do not count if ode has limited time
    else
        result.error = min(abs(distance_error));       % [m]
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
    
    result.t = t;
end

function result = swing_by(t_voo, initial_pos, initial_velocity_guess, venus_initial_pos, impulse_fraction)    
    V_exit_earth = initial_velocity_guess;            % [km/s]
    
    if nargin <= 3
        loadData;
        theta_venus_f = 3.1405;   % From 2B optim
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
    
    result.t = t(1:i-1);
end

function delta_v_swing_by = perform_swing_by(theta_entrada)
    loadData;
    t_voo = 0.1;
    v_entra_swing_by = [-5.1471      -37.727            0]; % [km/s]

    % Simula o swing by
    venus_initial_pos = R_venus_sun * [1,0,0] * Rz(3.1405);
    initial_pos = venus_initial_pos + SOI_venus * [1,0,0] * Rz(theta_entrada);
    venus_swing_by_pre = swing_by(...
        t_voo, ... % Flight time Mars-Venus [ut]
        initial_pos, ...    % Departure from Venus [m]
        v_entra_swing_by,... % Velocity to reach Mars orbit from Venus [km/s]
        venus_initial_pos... % Venus position before swing-by
    );
    venus_swing_by = swing_by(...
        venus_swing_by_pre.t(end), ... % Flight time Mars-Venus [ut]
        initial_pos, ...    % Departure from Venus [m]
        v_entra_swing_by,... % Velocity to reach Mars orbit from Venus [km/s]
        venus_initial_pos... % Venus position before swing-by
    );

    v_sai_swing_by = [0.51576      -40.787            0];
    
    if t_voo <= 0
        delta_v_swing_by = Inf;
    else
        delta_v_swing_by = norm(v_sai_swing_by - venus_swing_by.v_t2/1000);
    end
end

function x = optimize_swing_by()
    global shouldPlot;
    shouldPlot = 0;
    fun = @(x) perform_swing_by(x(1));
    x0 = [pi/2];
    options = optimset('MaxIter',1e4, 'MaxFunEvals', 1e4, 'TolFun', 1e-8, 'TolX', 1e-8);
    x = fminsearch(fun,x0, options);
    disp(x);
    disp(fun(x));
    shouldPlot = 1;
end