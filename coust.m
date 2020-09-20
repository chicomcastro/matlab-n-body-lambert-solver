function value = coust(x)
    t_voo = x(1)
    V_earth_mars = x(2)
    main;
    scPos = y(:,7:9);
    scVel = y(:,16:18); 
    scPosRadial = vecnorm(scPos')';
    distanceToMars = scPosRadial - R_mars_sun/ud;
    
    value = min(abs(distanceToMars));
end