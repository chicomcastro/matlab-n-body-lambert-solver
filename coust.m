function [posValue, vecValue] = coust(x)
    t_voo = x(1);           % ut
    V_earth_mars = x(2);    % km/s
    main;
    scPos = y(:,7:9);
    scVel = y(:,16:18); 
    scPosRadial = vecnorm(scPos')';
    distanceToMars = scPosRadial - R_mars_sun/ud;
    
    h = cross(scPos, scVel);
    v_versor = cross(h, scPos)./vecnorm(cross(h, scPos)')';
    V_mars_sun = sqrt(G_metric*M_sun/R_mars_sun)*v_versor/ud*ut;
    
    posValue = abs(distanceToMars(end, :)); % UA
    vecValue = norm(scVel(end, :) - V_mars_sun(end, :))*ud/ut/1000; % km/s
    
    %value = posValue + vecValue;
end