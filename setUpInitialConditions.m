%% Set up initial conditions

if ~exist("N", "var")
    setUpParameters;
end

%% Earth-Sun-Spacechip system
% Setting first body as the Sun
Mass(1) = M_sun/um;
r1 = zeros(1,3);
v1 = zeros(1,3);

Mass(2) = M_venus/um;
if ~exist("venus_initial_pos", "var")
    theta_venus_f = 2.14795059794092;   % From 2B optim
    omega_venus_sun = norm([0,0,3.236706097374289e-07])*ut; % [rad/ut]
    theta_venus_i = theta_venus_f - omega_venus_sun * t_voo;
    venus_initial_pos = R_venus_sun * [1,0,0] * Rz(theta_venus_i);
end
r2 = venus_initial_pos/ud;  % [ud]
tangent_direction = cross([0,0,1], r2(:)')/norm(cross([0,0,1], r2(:)'));
v2 = V_venus_sun/ud*ut * tangent_direction;  % [ud/ut]

% Spaceship (green)
Mass(3) = spaceship_mass/um;
if ~exist("initial_pos", "var")
    initial_pos = (R_earth_sun * [1 0 0] + 300e3 * [1 0 0]);  % m
end
if ~exist("V_exit_earth", "var")
    V_exit_earth = V_hohmann_earth_mars; %km/s
end
r3 = initial_pos/ud;
v3 = (V_exit_earth*1000)/ud*ut;

%%
% Sets initial conditions vector
y0 = [r1, r2, r3, v1, v2, v3];
y0 = y0(:);     % [column vector]