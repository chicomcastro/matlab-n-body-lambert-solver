function result = coust(x)
    t_voo = x.t_voo;           % ut
    initial_pos = x.initial_pos;
    V_exit_earth = x.v0;    % km/s
    
    main;
    
    r = y(end,7:9);           % ud
    v = y(end,16:18);         % ud/ut
    
    target = x.target_pos;        % ud
    error = r(:)' - target(:)';
    
    result.error = error;
    result.x_t1 = y(1,7:9)*ud;
    result.v_t1 = y(1,16:18)*ud/ut;
    result.x_t2 = r*ud;
    result.v_t2 = v*ud/ut;
end