
% Velocity spacial determination from r v    
h = cross(final_r, final_v);
v_versor = cross(h, final_r)./vecnorm(cross(h, final_r)')';
V_mars_sun = sqrt(G_metric*M_sun/R_mars_sun)*v_versor/ud*ut;

v = norm(scVel(end, :) - V_mars_sun(end, :))*ud/ut/1000; % km/s