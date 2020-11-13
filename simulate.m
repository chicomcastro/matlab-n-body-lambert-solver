function result = simulate(x)
% Run simulation for parameters passed on x argument and returns results
% on result struct. Assumes the error on distance as error for the solvers
    t_voo = x.t_voo;                % [ut]
    initial_pos = x.initial_pos;    % [m]
    V_exit_earth = x.v0;            % [km/s]
    
    if ~isfield(x, 'venus_initial_pos')
        loadData;
        theta_venus_f = 2.14795059794092;   % From 2B optim
        omega_venus_sun = norm([0,0,3.236706097374289e-07])*ut; % [rad/ut]
        theta_venus_i = theta_venus_f - omega_venus_sun * t_voo;
        venus_initial_pos = R_venus_sun * [1,0,0] * Rz(theta_venus_i);
    else
        venus_initial_pos = x.venus_initial_pos;
    end
    
    % Run all subroutines related to simulation
    main;
    
    global N;
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