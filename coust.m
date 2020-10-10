function result = coust(x)
% Run simulation for parameters passed on x argument and returns results
% on result struct. Assumes the error on distance as error for the solvers
    t_voo = x.t_voo;           % ut
    initial_pos = x.initial_pos;
    V_exit_earth = x.v0;    % km/s
    
    % Run all subroutines related to simulation
    main;
    
    r = y(end,7:9);                 % [ud]
    v = y(end,16:18);               % [ud/ut]
    
    target = x.target_pos;          % [ud]
    error = r(:)'*ud - target(:)';  % [m] [row vector]
    
    result.error = error;           % [m]
    result.x_t1 = y(1,7:9)*ud;      % [m]
    result.v_t1 = y(1,16:18)*ud/ut; % [m/s]
    result.x_t2 = r*ud;             % [m]
    result.v_t2 = v*ud/ut;          % [m/s]
end